# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests
import shutil

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from kb_orthofinder.kb_orthofinderImpl import kb_orthofinder
from kb_orthofinder.kb_orthofinderServer import MethodContext
from kb_orthofinder.authclient import KBaseAuth as _KBaseAuth

class kb_orthofinderTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_orthofinder'):
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
                            {'service': 'kb_orthofinder',
                             'method': 'annotate_plant_transcripts',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_orthofinder(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.test_data = cls.cfg['test_data']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.gfu = GenomeFileUtil(cls.callback_url)
        cls.dfu = DataFileUtil(cls.callback_url)
        cls.genome = "Test_Genome"
        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
#            print('Test workspace '+cls.wsName+' was not deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_orthofinder_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    @classmethod
    def prepare_data(cls):
        cls.gff_filename = 'Test_v1.0.gene.gff3.gz'
        cls.gff_path = os.path.join(cls.scratch, cls.gff_filename)
        shutil.copy(os.path.join("/kb", "module", "data", cls.gff_filename), cls.gff_path)

        cls.fa_filename = 'Test_v1.0.fa.gz'
        cls.fa_path = os.path.join(cls.scratch, cls.fa_filename)
        shutil.copy(os.path.join("/kb", "module", "data", cls.fa_filename), cls.fa_path)

        cls.tr_filename = cls.test_data+'.tar.gz'
        cls.tr_path = os.path.join("/kb", "module", "data", cls.tr_filename)
        # cls.tr_path = os.path.join(cls.scratch, cls.tr_filename)
        # This is now copied and unarchived in scripts/entrypoint.sh
        # shutil.copy(os.path.join("/kb", "module", "data", cls.tr_filename), cls.tr_path)

    def loadFakeGenome(cls):
        
        input_params = {
            'fasta_file': {'path': cls.fa_path},
            'gff_file': {'path': cls.gff_path},
            'genome_name': cls.genome,
            'workspace_name': cls.getWsName(),
            'source': 'Phytozome',
            'type': 'Reference',
            'scientific_name': 'Populus trichocarpa'
        }

        result = cls.gfu.fasta_gff_to_genome(input_params)

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_annotate_plant_transcripts(self):

        # Load Fake Genome
        self.loadFakeGenome()

        # DFU hangs on the large archive so we changed it so it's done
        # during test initialization in scripts/entrypoint.sh
        # self.dfu.unpack_file({'file_path' : self.tr_path})
        unpacked_tr = self.tr_path.replace('.tar.gz','')

        # Running Plant RAST
        ret = self.getImpl().annotate_plant_transcripts(self.getContext(), {'input_ws' : self.getWsName(),
                                                                            'input_genome' : self.genome,
                                                                            'families_path' : unpacked_tr,
                                                                            'threshold' : 0.55})

        print("RESULT: ",ret[0])
        self.assertEqual(ret[0]['transcripts'],2030)
        self.assertEqual(ret[0]['alignments'],1754)
        self.assertEqual(ret[0]['ftrs'],975)
        self.assertEqual(ret[0]['fns'],816)
        self.assertEqual(ret[0]['hit_ftrs'],73)
        self.assertEqual(ret[0]['hit_fns'],39)
        pass
