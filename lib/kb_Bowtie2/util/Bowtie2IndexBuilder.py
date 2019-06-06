from pprint import pprint

import os
import time
import traceback

from kb_Bowtie2.util.Bowtie2Runner import Bowtie2Runner

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.GenomeAnnotationAPIServiceClient import GenomeAnnotationAPI
from installed_clients.DataFileUtilClient import DataFileUtil


class Bowtie2IndexBuilder(object):

    def __init__(self, scratch_dir, ws_url, callback_url, service_wizard_url, provenance):
        self.scratch_dir = scratch_dir
        self.ws_url = ws_url
        self.ws = Workspace(self.ws_url)
        self.callback_url = callback_url
        self.service_wizard_url = service_wizard_url
        self.bowtie2 = Bowtie2Runner(self.scratch_dir)
        self.provenance = provenance


    def get_index(self, params):
        ''' The key function of this module- get a bowtie2 index for the specified input '''

        # validate the parameters and fetch assembly_info
        validated_params = self._validate_params(params)
        assembly_info = self._get_assembly_info(validated_params['ref'])

        # check the cache (keyed off of assembly_info)
        index_info = self._get_cached_index(assembly_info, validated_params)
        if index_info:
            index_info['from_cache'] = 1
            index_info['pushed_to_cache'] = 0
        else:
            # on a cache miss, build the index
            index_info = self._build_index(assembly_info, validated_params)
            index_info['from_cache'] = 0
            # pushed_to_cache will be set in return from _build_index

        index_info['assembly_ref'] = assembly_info['ref']
        index_info['genome_ref'] = assembly_info['genome_ref']

        return index_info


    def _validate_params(self, params):
        ''' validate parameters; can do some processing here to produce validated params '''
        validated_params = {'ref': None}
        if 'ref' in params and params['ref']:
            validated_params['ref'] = params['ref']
        else:
            raise ValueError('"ref" field indicating either an assembly or genome is required.')

        if 'output_dir' in params:
            validated_params['output_dir'] = params['output_dir']
        else:
            validated_params['output_dir'] = os.path.join(self.scratch_dir,
                                                          'bowtie2_index_' + str(int(time.time() * 100)))

        if os.path.exists(validated_params['output_dir']):
            raise('Output directory name specified (' + validated_params['output_dir'] +
                  ') already exists. Will not overwrite, so aborting.')

        if 'ws_for_cache' in params and params['ws_for_cache']:
            validated_params['ws_for_cache'] = params['ws_for_cache']
        else:
            print('WARNING: bowtie2 index if created will not be cached because "ws_for_cache" field not set')
            validated_params['ws_for_cache'] = None

        return validated_params


    def _get_assembly_info(self, ref):
        ''' given a ref to an assembly or genome, figure out the assembly and return its info '''
        info = self.ws.get_object_info3({'objects': [{'ref': ref}]})['infos'][0]
        obj_type = info[2]
        if obj_type.startswith('KBaseGenomeAnnotations.Assembly') or obj_type.startswith('KBaseGenomes.ContigSet'):
            return {'info': info, 'ref': ref, 'genome_ref': None}

        if obj_type.startswith('KBaseGenomes.Genome'):
            # we need to get the assembly for this genome
            ga = GenomeAnnotationAPI(self.service_wizard_url)
            assembly_ref = ga.get_assembly({'ref': ref})
            # using the path ensures we can access the assembly even if we don't have direct access
            ref_path = ref + ';' + assembly_ref
            info = self.ws.get_object_info3({'objects': [{'ref': ref_path}]})['infos'][0]
            return {'info': info, 'ref': ref_path, 'genome_ref': ref}

        raise ValueError('Input object was not of type: Assembly, ContigSet or Genome.  Cannot get Bowtie2 Index.')



    def _get_cached_index(self, assembly_info, validated_params):

        try:
            # note: list_reference_objects does not yet support reference paths, so we need to call
            # with the direct reference.  So we won't get a cache hit if you don't have direct access
            # to the assembly object right now (although you can still always build the assembly object)
            # Once this call supports paths, this should be changed to set ref = assembly_info['ref']
            info = assembly_info['info']
            ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
            objs = self.ws.list_referencing_objects([{'ref': ref}])[0]

            # iterate through each of the objects that reference the assembly
            bowtie2_indexes = []
            for o in objs:
                if o[2].startswith('KBaseRNASeq.Bowtie2IndexV2'):
                    bowtie2_indexes.append(o)

            # Nothing refs this assembly, so cache miss
            if len(bowtie2_indexes) == 0:
                return False

            # if there is more than one hit, get the most recent one
            # (obj_info[3] is the save_date timestamp (eg 2017-05-30T22:56:49+0000), so we can sort on that)
            bowtie2_indexes.sort(key=lambda x: x[3])
            bowtie2_index_info = bowtie2_indexes[-1]
            index_ref = str(bowtie2_index_info[6]) + '/' + str(bowtie2_index_info[0]) + '/' + str(bowtie2_index_info[4])

            # get the object data
            index_obj_data = self.ws.get_objects2({'objects': [{'ref': index_ref}]})['data'][0]['data']

            # download the handle object
            os.makedirs(validated_params['output_dir'])

            dfu = DataFileUtil(self.callback_url)
            dfu.shock_to_file({'file_path': os.path.join(validated_params['output_dir'], 'bt2_index.tar.gz'),
                               'handle_id': index_obj_data['handle']['hid'],
                               'unpack': 'unpack'})
            print('Cache hit: ')
            pprint(index_obj_data)
            return {'output_dir': validated_params['output_dir'],
                    'index_files_basename': index_obj_data['index_files_basename']}


        except Exception:
            # if we fail in saving the cached object, don't worry
            print('WARNING: exception encountered when trying to lookup in cache:')
            print(traceback.format_exc())
            print('END WARNING: exception encountered when trying to lookup in cache.')


        return None

    def _put_cached_index(self, assembly_info, index_files_basename, output_dir, ws_for_cache):

        if not ws_for_cache:
            print('WARNING: bowtie2 index cannot be cached because "ws_for_cache" field not set')
            return False

        try:
            dfu = DataFileUtil(self.callback_url)
            result = dfu.file_to_shock({'file_path': output_dir,
                                        'make_handle': 1,
                                        'pack': 'targz'})

            bowtie2_index = {'handle': result['handle'], 'size': result['size'],
                             'assembly_ref': assembly_info['ref'],
                             'index_files_basename': index_files_basename}

            ws = Workspace(self.ws_url)
            save_params = {'objects': [{'hidden': 1,
                                        'provenance': self.provenance,
                                        'name': os.path.basename(output_dir),
                                        'data': bowtie2_index,
                                        'type': 'KBaseRNASeq.Bowtie2IndexV2'
                                        }]
                           }
            if ws_for_cache.strip().isdigit():
                save_params['id'] = int(ws_for_cache)
            else:
                save_params['workspace'] = ws_for_cache.strip()
            save_result = ws.save_objects(save_params)
            print('Bowtie2IndexV2 cached to: ')
            pprint(save_result[0])
            return True

        except Exception:
            # if we fail in saving the cached object, don't worry
            print('WARNING: exception encountered when trying to cache the index files:')
            print(traceback.format_exc())
            print('END WARNING: exception encountered when trying to cache the index files')

        return False


    def _build_index(self, assembly_info, validated_params):
        # get the assembly as a fasta file using AssemblyUtil
        au = AssemblyUtil(self.callback_url)
        fasta_info = au.get_assembly_as_fasta({'ref': assembly_info['ref']})

        # make the target destination folder (check again it wasn't created yet)
        if os.path.exists(validated_params['output_dir']):
            raise('Output directory name specified (' + validated_params['output_dir'] +
                  ') already exists. Will not overwrite, so aborting.')
        os.makedirs(validated_params['output_dir'])

        # configure the command line args and run it
        cli_params = self._build_cli_params(fasta_info['path'], fasta_info['assembly_name'], validated_params)
        self.bowtie2.run('bowtie2-build', cli_params)
        index_info = {'output_dir': validated_params['output_dir'],
                      'index_files_basename': fasta_info['assembly_name']}

        # cache the result, mark if it worked or not
        cache_success = self._put_cached_index(assembly_info,
                                               fasta_info['assembly_name'],
                                               validated_params['output_dir'],
                                               validated_params['ws_for_cache'])
        if cache_success:
            index_info['pushed_to_cache'] = 1
        else:
            index_info['pushed_to_cache'] = 0

        return index_info


    def _build_cli_params(self, fasta_file_path, index_files_basename, validated_params):
        cli_params = []

        # always run in quiet mode
        cli_params.append('--quiet')

        # positional args: first the fasta path, then the base name used for the index files
        cli_params.append(fasta_file_path)
        cli_params.append(os.path.join(validated_params['output_dir'], index_files_basename))

        return cli_params
