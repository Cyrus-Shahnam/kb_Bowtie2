# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from installed_clients.WorkspaceClient import Workspace as workspaceService
from kb_Bowtie2.kb_Bowtie2Impl import kb_Bowtie2
from kb_Bowtie2.kb_Bowtie2Server import MethodContext
from kb_Bowtie2.authclient import KBaseAuth as _KBaseAuth

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.SetAPIServiceClient import SetAPI
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil


class kb_Bowtie2AlignerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_Bowtie2'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_Bowtie2',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_Bowtie2(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.srv_wiz_url = cls.cfg['srv-wiz-url']

        cls.dfu = DataFileUtil(cls.callback_url)
        cls.gfu = GenomeFileUtil(cls.callback_url)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        # return 'test_kb_Bowtie2_1499812859552'
        suffix = int(time.time() * 1000)
        wsName = "test_kb_Bowtie2_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def loadSingleEndReads(self):
        if hasattr(self.__class__, 'se_reads_ref'):
            return self.__class__.se_reads_ref
        # return '23735/2/1'
        fq_path = os.path.join(self.scratch, 'extracted_WT_rep1.fastq')
        shutil.copy(os.path.join('data', 'bt_test_data', 'extracted_WT_rep1.fastq'), fq_path)

        ru = ReadsUtils(self.callback_url)
        se_reads_ref = ru.upload_reads({'fwd_file': fq_path,
                                        'wsname': self.getWsName(),
                                        'name': 'test_readsSE',
                                        'sequencing_tech': 'artificial reads'})['obj_ref']
        self.__class__.se_reads_ref = se_reads_ref
        print('Loaded SingleEndReads: ' + se_reads_ref)
        return se_reads_ref

    def loadPairedEndReads(self):
        if hasattr(self.__class__, 'pe_reads_ref'):
            return self.__class__.pe_reads_ref
        # return '23735/3/1'
        fq_path1 = os.path.join(self.scratch, 'reads_1.fq')
        shutil.copy(os.path.join('data', 'bt_test_data', 'reads_1.fq'), fq_path1)
        fq_path2 = os.path.join(self.scratch, 'reads_2.fq')
        shutil.copy(os.path.join('data', 'bt_test_data', 'reads_2.fq'), fq_path2)

        ru = ReadsUtils(self.callback_url)
        pe_reads_ref = ru.upload_reads({'fwd_file': fq_path1, 'rev_file': fq_path2,
                                        'wsname': self.getWsName(),
                                        'name': 'test_readsPE',
                                        'sequencing_tech': 'artificial reads'})['obj_ref']
        self.__class__.pe_reads_ref = pe_reads_ref
        print('Loaded PairedEndReads: ' + pe_reads_ref)
        return pe_reads_ref

    def loadGenome(self):
        if hasattr(self.__class__, 'genome_ref'):
            return self.__class__.genome_ref

        test_gbk_file = os.path.join("data", "bt_test_data", "at_chrom1_section.gbk")
        gbk_file = os.path.join(self.scratch, os.path.basename(test_gbk_file))
        shutil.copy(test_gbk_file, gbk_file)

        genome_ref = self.gfu.genbank_to_genome({
            "file": {
                "path": gbk_file
            },
            "genome_name": 'my_test_genome',
            "workspace_name": self.getWsName(),
            "source": "Ensembl",
            "type": "User upload",
            "generate_ids_if_needed": 1
        }).get('genome_ref')

        self.__class__.genome_ref = genome_ref
        print('Loaded Genome: ' + genome_ref)
        return genome_ref

    def loadAssembly(self):
        if hasattr(self.__class__, 'assembly_ref'):
            return self.__class__.assembly_ref
        # return '23735/1/1'
        fasta_path = os.path.join(self.scratch, 'test_ref.fa')
        shutil.copy(os.path.join('data', 'bt_test_data', 'test_ref.fa'), fasta_path)
        au = AssemblyUtil(self.callback_url)
        assembly_ref = au.save_assembly_from_fasta({'file': {'path': fasta_path},
                                                    'workspace_name': self.getWsName(),
                                                    'assembly_name': 'test_assembly'
                                                    })
        self.__class__.assembly_ref = assembly_ref
        print('Loaded Assembly: ' + assembly_ref)
        return assembly_ref

    def loadSampleSet(self):
        if hasattr(self.__class__, 'sample_set_ref'):
            return self.__class__.sample_set_ref
        # return '23735/4/1'
        pe_reads_ref = self.loadPairedEndReads()
        sample_set_name = 'TestSampleSet'
        sample_set_data = {'Library_type': 'PairedEnd',
                           'domain': "Prokaryotes",
                           'num_samples': 3,
                           'platform': None,
                           'publication_id': None,
                           'sample_ids': [pe_reads_ref, pe_reads_ref, pe_reads_ref],
                           'sampleset_desc': None,
                           'sampleset_id': sample_set_name,
                           'condition': ['ss1', 'ss2', 'ss3'],
                           'source': None
                           }
        ss_obj = self.getWsClient().save_objects({'workspace': self.getWsName(),
                                                  'objects': [{'type': 'KBaseRNASeq.RNASeqSampleSet',
                                                               'data': sample_set_data,
                                                               'name': sample_set_name,
                                                               'provenance': [{}]
                                                               }]
                                                  })
        ss_ref = "{}/{}/{}".format(ss_obj[0][6], ss_obj[0][0], ss_obj[0][4])
        print('Loaded SampleSet: ' + ss_ref)
        return ss_ref

    def loadReadsSet(self):
        if hasattr(self.__class__, 'reads_set_ref'):
            return self.__class__.reads_set_ref
        pe_reads_ref = self.loadPairedEndReads()
        reads_set_name = 'TestReadsSet'
        # create the set object

        reads_set_data = {
            'description': 'Reads Set for testing Bowtie',
            'items': [{
                'ref': pe_reads_ref,
                'label': 'rs1'
            }, {
                'ref': pe_reads_ref,
                'label': 'rs2'
            }, {
                'ref': pe_reads_ref,
                'label': 'rs3'
            }
            ]
        }
        # test a save
        set_api = SetAPI(self.srv_wiz_url)
        res = set_api.save_reads_set_v1({
            'data': reads_set_data,
            'output_object_name': reads_set_name,
            'workspace': self.getWsName()
        })
        reads_set_ref = res['set_ref']

        # reads_set_ref = '5264/52/1'
        print('Loaded ReadsSet: ' + reads_set_ref)
        return reads_set_ref

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_bowtie2_aligner_with_sampleset(self):
        assembly_ref = self.loadAssembly()
        genome_ref = self.loadGenome()
        se_lib_ref = self.loadSingleEndReads()
        params = {'input_ref': se_lib_ref,
                  'assembly_or_genome_ref': genome_ref,
                  'output_obj_name_suffix': 'readsAlignment1',
                  'output_alignment_suffix': '_some_ext',
                  'output_workspace': self.getWsName(),
                  'concurrent_njsw_tasks': 0,
                  'concurrent_local_tasks': 1}
        pprint(params)
        res = self.getImpl().align_reads_to_assembly_app(self.getContext(), params)[0]
        pprint(res)
        self.assertIn('report_info', res)
        self.assertIn('report_name', res['report_info'])
        self.assertIn('report_ref', res['report_info'])

        obj_ref = res.get('output_info').get('upload_results').get('obj_ref')
        obj_data = self.dfu.get_objects(
            {"object_refs": [obj_ref]})['data'][0]['data']
        align_stats = obj_data.get('alignment_stats')
        self.assertEqual(align_stats.get('total_reads'), 15254)
        self.assertEqual(align_stats.get('mapped_reads'), 15205)
        self.assertEqual(align_stats.get('unmapped_reads'), 49)
        self.assertEqual(align_stats.get('singletons'), 4835)
        self.assertEqual(align_stats.get('multiple_alignments'), 10370)

        ss_ref = self.loadSampleSet()
        params = {'input_ref': ss_ref,
                  'assembly_or_genome_ref': assembly_ref,
                  'output_obj_name_suffix': 'readsAlignment1',
                  'output_alignment_suffix': '_some_ext',
                  'output_workspace': self.getWsName(),
                  'concurrent_njsw_tasks': 0,
                  'concurrent_local_tasks': 1}
        pprint('Running with a SampleSet')
        pprint(params)
        res = self.getImpl().align_reads_to_assembly_app(self.getContext(), params)[0]
        pprint(res)
        self.assertIn('report_info', res)
        self.assertIn('report_name', res['report_info'])
        self.assertIn('report_ref', res['report_info'])

    def test_bowtie2_aligner_with_readsset(self):
        assembly_ref = self.loadAssembly()
        se_lib_ref = self.loadSingleEndReads()
        params = {'input_ref': se_lib_ref,
                  'assembly_or_genome_ref': assembly_ref,
                  'output_obj_name_suffix': 'readsAlignment1',
                  'output_alignment_suffix': '_some_ext',
                  'output_workspace': self.getWsName(),
                  'concurrent_njsw_tasks': 0,
                  'concurrent_local_tasks': 1}
        pprint(params)
        res = self.getImpl().align_reads_to_assembly_app(self.getContext(), params)[0]
        pprint(res)
        self.assertIn('report_info', res)
        self.assertIn('report_name', res['report_info'])
        self.assertIn('report_ref', res['report_info'])

        reads_set_ref = self.loadReadsSet()
        params = {'input_ref': reads_set_ref,
                  'assembly_or_genome_ref': assembly_ref,
                  'output_obj_name_suffix': 'readsAlignment1',
                  'output_alignment_suffix': '_some_ext',
                  'output_workspace': self.getWsName(),
                  'concurrent_njsw_tasks': 0,
                  'concurrent_local_tasks': 1}
        pprint('Running with a ReadsSet')
        pprint(params)
        res = self.getImpl().align_reads_to_assembly_app(self.getContext(), params)[0]
        pprint(res)
        self.assertIn('report_info', res)
        self.assertIn('report_name', res['report_info'])
        self.assertIn('report_ref', res['report_info'])
