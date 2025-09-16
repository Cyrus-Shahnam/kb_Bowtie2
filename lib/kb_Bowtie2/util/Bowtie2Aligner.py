from pprint import pprint

import os
import time
import uuid
import copy
import re

from kb_Bowtie2.util.Bowtie2Runner import Bowtie2Runner
from kb_Bowtie2.util.Bowtie2IndexBuilder import Bowtie2IndexBuilder

from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from installed_clients.kb_QualiMapClient import kb_QualiMap
from installed_clients.KBaseReportClient import KBaseReport

from installed_clients.WorkspaceClient import Workspace
from installed_clients.SetAPIServiceClient import SetAPI

from installed_clients.KBParallelClient import KBParallel


class Bowtie2Aligner(object):

    def __init__(self, scratch_dir, workspace_url, callback_url, srv_wiz_url, provenance):
        self.scratch_dir = scratch_dir
        self.workspace_url = workspace_url
        self.callback_url = callback_url
        self.srv_wiz_url = srv_wiz_url
        self.provenance = provenance

        # from the provenance, extract out the version to run by exact hash if possible
        self.my_version = 'release'
        if len(provenance) > 0:
            if 'subactions' in provenance[0]:
                self.my_version = self.get_version_from_subactions('kb_Bowtie2', provenance[0]['subactions'])
        print('Running kb_Bowtie2 version = ' + self.my_version)

        self.ws = Workspace(self.workspace_url)
        self.bowtie2 = Bowtie2Runner(self.scratch_dir)
        self.parallel_runner = KBParallel(self.callback_url)
        self.qualimap = kb_QualiMap(self.callback_url)

    def get_version_from_subactions(self, module_name, subactions):
        # go through each sub action looking for
        if not subactions:
            return 'release'  # default to release if we can't find anything
        for sa in subactions:
            if 'name' in sa:
                if sa['name'] == module_name:
                    # local-docker-image implies that we are running in kb-test, so return 'dev'
                    if sa['commit'] == 'local-docker-image':
                        return 'dev'
                    # to check that it is a valid hash, make sure it is the right
                    # length and made up of valid hash characters
                    if re.match('[a-fA-F0-9]{40}$', sa['commit']):
                        return sa['commit']
        # again, default to setting this to release
        return 'release'

    def align(self, params):
        validated_params = self.validate_params(params)
        input_info = self.determine_input_info(validated_params)
        # input info provides information on the input and tells us if we should
        # run as a single_library or as a set:
        #     input_info = {'run_mode': '', 'info': [..], 'ref': '55/1/2'}

        assembly_or_genome_ref = validated_params['assembly_or_genome_ref']

        if input_info['run_mode'] == 'single_library':
            if 'output_alignment_name' not in validated_params:
                suffix = '_alignment'
                if 'output_alignment_suffix' in validated_params:
                    suffix = validated_params['output_alignment_suffix']
                validated_params['output_alignment_name'] = input_info['info'][1] + suffix
            single_lib_result = self.single_reads_lib_run(
                input_info,
                assembly_or_genome_ref,
                validated_params,
                create_report=validated_params['create_report']
            )

            return single_lib_result

        if input_info['run_mode'] == 'sample_set':
            reads = self.fetch_reads_refs_from_sampleset(input_info['ref'], input_info['info'], validated_params)
            # pass use_sais through so the index builder can add --sais if requested
            self.build_bowtie2_index(assembly_or_genome_ref, validated_params['output_workspace'], validated_params)

            print('Running on set of reads=')
            pprint(reads)

            tasks = []
            for r in reads:
                tasks.append(self.build_single_execution_task(r['ref'], params, r['alignment_output_name'], r['condition']))

            batch_run_params = {'tasks': tasks,
                                'runner': 'parallel',
                                'max_retries': 2}
            if validated_params['concurrent_local_tasks'] is not None:
                batch_run_params['concurrent_local_tasks'] = validated_params['concurrent_local_tasks']
            if validated_params['concurrent_njsw_tasks'] is not None:
                batch_run_params['concurrent_njsw_tasks'] = validated_params['concurrent_njsw_tasks']
            results = self.parallel_runner.run_batch(batch_run_params)
            print('Batch run results=')
            pprint(results)
            batch_result = self.process_batch_result(results, validated_params, reads, input_info['info'])
            return batch_result

        raise ('Improper run mode')

    def build_single_execution_task(self, reads_lib_ref, params, output_name, condition):
        task_params = copy.deepcopy(params)

        task_params['input_ref'] = reads_lib_ref
        task_params['output_alignment_name'] = output_name
        task_params['create_report'] = 0
        task_params['condition_label'] = condition

        return {'module_name': 'kb_Bowtie2',
                'function_name': 'align_reads_to_assembly_app',
                'version': self.my_version,
                'parameters': task_params}

    def single_reads_lib_run(self, read_lib_info, assembly_or_genome_ref, validated_params,
                             create_report=False, bowtie2_index_info=None):
        ''' run on one reads '''

        # download reads and prepare any bowtie2 index files
        input_configuration = self.prepare_single_run(
            read_lib_info, assembly_or_genome_ref,
            bowtie2_index_info, validated_params['output_workspace'],
            validated_params
        )

        # run the actual program
        run_output_info = self.run_bowtie2_align_cli(input_configuration, validated_params)

        # process the result and save the output
        upload_results = self.save_read_alignment_output(run_output_info, input_configuration, validated_params)
        run_output_info['upload_results'] = upload_results

        report_info = None
        if create_report:
            report_info = self.create_report_for_single_run(run_output_info, input_configuration, validated_params)

        self.clean(run_output_info)

        return {'output_info': run_output_info, 'report_info': report_info}

    def build_bowtie2_index(self, assembly_or_genome_ref, ws_for_cache, validated_params=None):
        bowtie2IndexBuilder = Bowtie2IndexBuilder(self.scratch_dir, self.workspace_url,
                                                  self.callback_url, self.srv_wiz_url,
                                                  self.provenance)

        index_params = {'ref': assembly_or_genome_ref,
                        'ws_for_cache': ws_for_cache}
        # propagate SAIS preference if present
        if validated_params and 'use_sais' in validated_params:
            index_params['use_sais'] = validated_params['use_sais']

        return bowtie2IndexBuilder.get_index(index_params)

    def prepare_single_run(self, input_info, assembly_or_genome_ref,
                           bowtie2_index_info, ws_for_cache, validated_params):
        ''' Given a reads ref and an assembly, setup the bowtie2 index '''
        # first setup the bowtie2 index of the assembly
        input_configuration = {'bowtie2_index_info': bowtie2_index_info}
        if not bowtie2_index_info:
            bowtie2IndexBuilder = Bowtie2IndexBuilder(self.scratch_dir, self.workspace_url,
                                                      self.callback_url, self.srv_wiz_url,
                                                      self.provenance)

            index_params = {'ref': assembly_or_genome_ref,
                            'ws_for_cache': ws_for_cache}
            # propagate SAIS preference if present
            if 'use_sais' in validated_params:
                index_params['use_sais'] = validated_params['use_sais']

            index_result = bowtie2IndexBuilder.get_index(index_params)
            input_configuration['bowtie2_index_info'] = index_result

        # next download the reads
        read_lib_ref = input_info['ref']
        read_lib_info = input_info['info']
        reads_params = {'read_libraries': [read_lib_ref],
                        'interleaved': 'false',
                        'gzipped': None}
        ru = ReadsUtils(self.callback_url)
        reads = ru.download_reads(reads_params)['files']

        input_configuration['reads_lib_type'] = self.get_type_from_obj_info(read_lib_info).split('.')[1]
        input_configuration['reads_files'] = reads[read_lib_ref]
        input_configuration['reads_lib_ref'] = read_lib_ref

        return input_configuration

    def run_bowtie2_align_cli(self, input_configuration, validated_params):
        # pprint('======== input_configuration =====')
        # pprint(input_configuration)
        options = []
        run_output_info = {}

        # set the bowtie2 index location
        bt2_index_dir = input_configuration['bowtie2_index_info']['output_dir']
        bt2_index_basename = input_configuration['bowtie2_index_info']['index_files_basename']
        options.extend(['-x', bt2_index_basename])

        # set the input reads
        if input_configuration['reads_lib_type'] == 'SingleEndLibrary':
            options.extend(['-U', input_configuration['reads_files']['files']['fwd']])
            run_output_info['library_type'] = 'single_end'
        elif input_configuration['reads_lib_type'] == 'PairedEndLibrary':
            options.extend(['-1', input_configuration['reads_files']['files']['fwd']])
            options.extend(['-2', input_configuration['reads_files']['files']['rev']])
            run_output_info['library_type'] = 'paired_end'

        # setup the output file name
        output_dir = os.path.join(self.scratch_dir, 'bowtie2_alignment_output_' + str(int(time.time() * 10000)))
        output_sam_file = os.path.join(output_dir, 'reads_alignment.sam')
        os.makedirs(output_dir)
        options.extend(['-S', output_sam_file])
        run_output_info['output_sam_file'] = output_sam_file
        run_output_info['output_dir'] = output_dir

        # parse all the other parameters
        if 'quality_score' in validated_params and str(validated_params['quality_score']).strip():
            options.append('--' + str(validated_params['quality_score']).strip())

        if 'alignment_type' in validated_params and str(validated_params['alignment_type']).strip():
            options.append('--' + str(validated_params['alignment_type']).strip())

        if 'preset_options' in validated_params:
            preset = str(validated_params['preset_options']).strip()
            if preset:
                if 'alignment_type' in validated_params and validated_params['alignment_type'] == 'local':
                    options.append('--' + preset + '-local')
                else:
                    options.append('--' + preset)

        if 'trim5' in validated_params:
            options.extend(['--trim5', str(validated_params['trim5'])])
        if 'trim3' in validated_params:
            options.extend(['--trim3', str(validated_params['trim3'])])
        if 'np' in validated_params:
            options.extend(['--np', str(validated_params['np'])])

        if 'minins' in validated_params:
            options.extend(['--minins', str(validated_params['minins'])])
        if 'maxins' in validated_params:
            options.extend(['--maxins', str(validated_params['maxins'])])

        # NEW: pass-through for Bowtie2 >= 2.5 sam-opt-config
        if 'sam_opt_config' in validated_params:
            sam_cfg = str(validated_params['sam_opt_config']).strip()
            if sam_cfg:
                options.extend(['--sam-opt-config', sam_cfg])

        options.extend(['-a'])  # report all distinct reads

        # unfortunately, bowtie2 expects the index files to be in the current directory, and
        # you cannot configure it otherwise.  So run bowtie out of the index directory, but
        # place the output SAM file somewhere else
        self.bowtie2.run('bowtie2', options, cwd=bt2_index_dir)

        return run_output_info

    def save_read_alignment_output(self, run_output_info, input_configuration, validated_params):
        rau = ReadsAlignmentUtils(self.callback_url)
        destination_ref = validated_params['output_workspace'] + '/' + validated_params['output_alignment_name']
        condition = 'unknown'
        if 'condition_label' in validated_params:
            condition = validated_params['condition_label']
        upload_params = {'file_path': run_output_info['output_sam_file'],
                         'destination_ref': destination_ref,
                         'read_library_ref': input_configuration['reads_lib_ref'],
                         'assembly_or_genome_ref': validated_params['assembly_or_genome_ref'],
                         'condition': condition}
        upload_results = rau.upload_alignment(upload_params)
        return upload_results

    def clean(self, run_output_info):
        ''' Not really necessary on a single run, but if we are running multiple local subjobs, we
        should clean up files that have already been saved back up to kbase '''
        pass

    def create_report_for_single_run(self, run_output_info, input_configuration, validated_params):
        # first run qualimap
        qualimap_report = self.qualimap.run_bamqc({'input_ref': run_output_info['upload_results']['obj_ref']})
        qc_result_zip_info = qualimap_report['qc_result_zip_info']

        # create report
        report_text = 'Ran on a single reads library.\n\n'
        alignment_info = self.get_obj_info(run_output_info['upload_results']['obj_ref'])
        report_text = 'Created ReadsAlignment: ' + str(alignment_info[1]) + '\n'
        report_text = '                        ' + run_output_info['upload_results']['obj_ref'] + '\n'
        kbr = KBaseReport(self.callback_url)
        report_info = kbr.create_extended_report({'message': report_text,
                                                  'objects_created': [{'ref': run_output_info['upload_results']['obj_ref'],
                                                                       'description': 'ReadsAlignment'}],
                                                  'report_object_name': 'kb_Bowtie2_' + str(uuid.uuid4()),
                                                  'direct_html_link_index': 0,
                                                  'html_links': [{'shock_id': qc_result_zip_info['shock_id'],
                                                                  'name': qc_result_zip_info['index_html_file_name'],
                                                                  'label': qc_result_zip_info['name']}],
                                                  'workspace_name': validated_params['output_workspace']
                                                  })
        return {'report_name': report_info['name'], 'report_ref': report_info['ref']}

    def process_batch_result(self, batch_result, validated_params, reads, input_set_info):

        n_jobs = len(batch_result['results'])
        n_success = 0
        n_error = 0
        ran_locally = 0
        ran_njsw = 0

        # reads alignment set items
        items = []
        objects_created = []

        for k in range(0, len(batch_result['results'])):
            job = batch_result['results'][k]
            result_package = job['result_package']
            if job['is_error']:
                n_error += 1
            else:
                n_success += 1
                output_info = result_package['result'][0]['output_info']
                ra_ref = output_info['upload_results']['obj_ref']
                # Note: could add a label to the alignment here?
                items.append({'ref': ra_ref, 'label': reads[k]['condition']})
                objects_created.append({'ref': ra_ref})

            if result_package['run_context']['location'] == 'local':
                ran_locally += 1
            if result_package['run_context']['location'] == 'njsw':
                ran_njsw += 1

        # Save the alignment set
        alignment_set_data = {'description': '', 'items': items}
        alignment_set_save_params = {'data': alignment_set_data,
                                     'workspace': validated_params['output_workspace'],
                                     'output_object_name': str(input_set_info[1]) + validated_params['output_obj_name_suffix']}

        set_api = SetAPI(self.srv_wiz_url)
        save_result = set_api.save_reads_alignment_set_v1(alignment_set_save_params)
        print('Saved ReadsAlignment=')
        pprint(save_result)
        objects_created.append({'ref': save_result['set_ref'], 'description': 'Set of all reads alignments generated'})
        set_name = save_result['set_info'][1]

        # run qualimap
        qualimap_report = self.qualimap.run_bamqc({'input_ref': save_result['set_ref']})
        qc_result_zip_info = qualimap_report['qc_result_zip_info']

        # create the report
        report_text = 'Ran on SampleSet or ReadsSet.\n\n'
        report_text = 'Created ReadsAlignmentSet: ' + str(set_name) + '\n\n'
        report_text += 'Total ReadsLibraries = ' + str(n_jobs) + '\n'
        report_text += '        Successful runs = ' + str(n_success) + '\n'
        report_text += '            Failed runs = ' + str(n_error) + '\n'
        report_text += '       Ran on main node = ' + str(ran_locally) + '\n'
        report_text += '   Ran on remote worker = ' + str(ran_njsw) + '\n\n'

        print('Report text=')
        print(report_text)

        kbr = KBaseReport(self.callback_url)
        report_info = kbr.create_extended_report({'message': report_text,
                                                  'objects_created': objects_created,
                                                  'report_object_name': 'kb_Bowtie2_' + str(uuid.uuid4()),
                                                  'direct_html_link_index': 0,
                                                  'html_links': [{'shock_id': qc_result_zip_info['shock_id'],
                                                                  'name': qc_result_zip_info['index_html_file_name'],
                                                                  'label': qc_result_zip_info['name']}],
                                                  'workspace_name': validated_params['output_workspace']
                                                  })

        result = {'report_info': {'report_name': report_info['name'], 'report_ref': report_info['ref']}}
        result['batch_output_info'] = batch_result

        return result

    def validate_params(self, params):
        validated_params = {}

        required_string_fields = ['input_ref', 'assembly_or_genome_ref', 'output_obj_name_suffix',
                                  'output_workspace']
        for field in required_string_fields:
            if field in params and params[field]:
                validated_params[field] = params[field]
            else:
                raise ValueError('"' + field + '" field required to run bowtie2 aligner app')

        optional_fields = [
            'quality_score', 'alignment_type', 'preset_options', 'trim5', 'trim3', 'condition_label',
            'np', 'minins', 'maxins', 'output_alignment_suffix', 'output_alignment_name',
            # NEW:
            'sam_opt_config', 'use_sais'
        ]
        for field in optional_fields:
            if field in params:
                if params[field] is not None:
                    # For strings, avoid passing empty values where they produce bad flags
                    if isinstance(params[field], str) and not params[field].strip():
                        if field in ('preset_options', 'sam_opt_config'):
                            continue
                    validated_params[field] = params[field]

        validated_params['create_report'] = True
        if 'create_report' in params and params['create_report'] is not None:
            if int(params['create_report']) == 1:
                validated_params['create_report'] = True
            elif int(params['create_report']) == 0:
                validated_params['create_report'] = False
            else:
                raise ValueError('"create_report" field, if present, should be set to a boolean value: 0 or 1')

        validated_params['concurrent_local_tasks'] = None
        validated_params['concurrent_njsw_tasks'] = None

        if 'concurrent_local_tasks' in params and params['concurrent_local_tasks'] is not None:
            validated_params['concurrent_local_tasks'] = int(params['concurrent_local_tasks'])
        if 'concurrent_njsw_tasks' in params and params['concurrent_njsw_tasks'] is not None:
            validated_params['concurrent_njsw_tasks'] = int(params['concurrent_njsw_tasks'])

        return validated_params

    def fetch_reads_refs_from_sampleset(self, ref, info, validated_params):
        """
        Note: adapted from kbaseapps/kb_hisat2 - file_util.py
        ...
        """
        obj_type = self.get_type_from_obj_info(info)
        refs = list()
        refs_for_ws_info = list()
        if "KBaseSets.ReadsSet" in obj_type or "KBaseRNASeq.RNASeqSampleSet" in obj_type:
            print("Looking up reads references in ReadsSet object")
            set_api = SetAPI(self.srv_wiz_url)
            reads_set = set_api.get_reads_set_v1({'ref': ref,
                                                  'include_item_info': 0,
                                                  'include_set_item_ref_paths': 1
                                                  })

            for reads in reads_set["data"]["items"]:
                refs.append({'ref': reads['ref_path'],
                             'condition': reads['label']
                             })
                refs_for_ws_info.append({'ref': reads['ref_path']})
        else:
            raise ValueError("Unable to fetch reads reference from object {} "
                             "which is a {}".format(ref, obj_type))

        # get object info so we can name things properly
        infos = self.ws.get_object_info3({'objects': refs_for_ws_info})['infos']

        name_ext = '_alignment'
        if 'output_alignment_suffix' in validated_params \
                and validated_params['output_alignment_suffix'] is not None:
            ext = validated_params['output_alignment_suffix'].replace(' ', '')
            if ext:
                name_ext = ext

        unique_name_lookup = {}
        for k in range(0, len(refs)):
            refs[k]['info'] = infos[k]
            name = infos[k][1]
            if name not in unique_name_lookup:
                unique_name_lookup[name] = 1
            else:
                unique_name_lookup[name] += 1
                name = name + '_' + str(unique_name_lookup[name])
            name = name + name_ext
            refs[k]['alignment_output_name'] = name

        return refs

    def determine_input_info(self, validated_params):
        ''' get info on the input_ref object and determine if we run once or run on a set '''
        info = self.get_obj_info(validated_params['input_ref'])
        obj_type = self.get_type_from_obj_info(info)
        if obj_type in ['KBaseAssembly.PairedEndLibrary', 'KBaseAssembly.SingleEndLibrary',
                        'KBaseFile.PairedEndLibrary', 'KBaseFile.SingleEndLibrary']:
            return {'run_mode': 'single_library', 'info': info, 'ref': validated_params['input_ref']}
        if obj_type == 'KBaseRNASeq.RNASeqSampleSet':
            return {'run_mode': 'sample_set', 'info': info, 'ref': validated_params['input_ref']}
        if obj_type == 'KBaseSets.ReadsSet':
            return {'run_mode': 'sample_set', 'info': info, 'ref': validated_params['input_ref']}

        raise ValueError('Object type of input_ref is not valid, was: ' + str(obj_type))

    def get_type_from_obj_info(self, info):
        return info[2].split('-')[0]

    def get_obj_info(self, ref):
        return self.ws.get_object_info3({'objects': [{'ref': ref}]})['infos'][0]
