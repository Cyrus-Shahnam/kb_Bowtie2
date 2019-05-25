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
from installed_clients.GenomeFileUtilClient import GenomeFileUtil


class kb_Bowtie2IndexerTest(unittest.TestCase):

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
        suffix = int(time.time() * 1000)
        wsName = "test_kb_Bowtie2_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName


    def loadAssembly(self):
        if hasattr(self.__class__, 'assembly_ref'):
            return self.__class__.assembly_ref
        fasta_path = os.path.join(self.scratch, 'test.fna')
        shutil.copy(os.path.join('data', 'test.fna'), fasta_path)
        au = AssemblyUtil(self.callback_url)
        assembly_ref = au.save_assembly_from_fasta({'file': {'path': fasta_path},
                                                    'workspace_name': self.getWsName(),
                                                    'assembly_name': 'test_assembly'
                                                    })
        self.__class__.assembly_ref = assembly_ref
        return assembly_ref


    def loadGenome(self):
        if hasattr(self.__class__, 'genome_ref'):
            return self.__class__.genome_ref
        genbank_file_path = os.path.join(self.scratch, 'minimal.gbff')
        shutil.copy(os.path.join('data', 'minimal.gbff'), genbank_file_path)
        gfu = GenomeFileUtil(self.callback_url)
        genome_ref = gfu.genbank_to_genome({'file': {'path': genbank_file_path},
                                            'workspace_name': self.getWsName(),
                                            'genome_name': 'test_genome',
                                            'source': 'Ensembl',
                                            'generate_ids_if_needed': 1,
                                            'generate_missing_genes': 1
                                            })['genome_ref']
        self.__class__.genome_ref = genome_ref
        return genome_ref


    def loadPEReads(self):
        if hasattr(self.__class__, 'assembly_ref'):
            return self.__class__.assembly_ref
        fasta_path = os.path.join(self.scratch, 'test.fna')
        shutil.copy(os.path.join('data', 'test.fna'), fasta_path)
        au = AssemblyUtil(self.callback_url)
        assembly_ref = au.save_assembly_from_fasta({'file': {'path': fasta_path},
                                                    'workspace_name': self.getWsName(),
                                                    'assembly_name': 'test_assembly'
                                                    })
        self.__class__.assembly_ref = assembly_ref
        return assembly_ref


    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    def test_bowtie2_installation(self):
        params = {'command': 'bowtie2', 'options': ['--help']}
        self.getImpl().run_bowtie2_cli(self.getContext(), params)


    def test_build_bowtie2_index_from_assembly(self):

        # test build directly from an assembly, forget to add ws_for_cache so object will not be cached
        assembly_ref = self.loadAssembly()
        res = self.getImpl().get_bowtie2_index(self.getContext(), {'ref': assembly_ref})[0]
        self.assertIn('output_dir', res)
        self.assertIn('from_cache', res)
        self.assertEquals(res['from_cache'], 0)
        self.assertIn('pushed_to_cache', res)
        self.assertEquals(res['pushed_to_cache'], 0)
        self.assertIn('index_files_basename', res)
        self.assertEquals(res['index_files_basename'], 'test_assembly')

        pprint(res)

        # do it again, and set ws_for_cache
        assembly_ref = self.loadAssembly()
        res = self.getImpl().get_bowtie2_index(self.getContext(), {'ref': assembly_ref,
                                                                   'ws_for_cache': self.getWsName()})[0]
        self.assertIn('output_dir', res)
        self.assertIn('from_cache', res)
        self.assertEquals(res['from_cache'], 0)
        self.assertIn('pushed_to_cache', res)
        self.assertEquals(res['pushed_to_cache'], 1)
        self.assertIn('index_files_basename', res)
        self.assertEquals(res['index_files_basename'], 'test_assembly')

        pprint(res)

        # do it again, should retrieve from cache
        assembly_ref = self.loadAssembly()
        res = self.getImpl().get_bowtie2_index(self.getContext(), {'ref': assembly_ref})[0]
        self.assertIn('output_dir', res)
        self.assertIn('from_cache', res)
        self.assertEquals(res['from_cache'], 1)
        self.assertIn('pushed_to_cache', res)
        self.assertEquals(res['pushed_to_cache'], 0)
        self.assertIn('index_files_basename', res)
        self.assertEquals(res['index_files_basename'], 'test_assembly')
        pprint(res)


    def test_build_bowtie2_index_from_genome(self):

        # finally, try it with a genome_ref instead
        genome_ref = self.loadGenome()
        res = self.getImpl().get_bowtie2_index(self.getContext(), {'ref': genome_ref})[0]
        self.assertIn('output_dir', res)
        self.assertIn('from_cache', res)
        self.assertEquals(res['from_cache'], 0)
        self.assertIn('pushed_to_cache', res)
        self.assertEquals(res['pushed_to_cache'], 0)
        self.assertIn('index_files_basename', res)
        self.assertEquals(res['index_files_basename'], 'test_genome_assembly')
        pprint(res)
