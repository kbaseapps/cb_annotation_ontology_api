# -*- coding: utf-8 -*-
import os
import time
import json
import unittest
from configparser import ConfigParser

from cb_annotation_ontology_api.cb_annotation_ontology_apiImpl import cb_annotation_ontology_api
from cb_annotation_ontology_api.cb_annotation_ontology_apiServer import MethodContext
from cb_annotation_ontology_api.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class cb_annotation_ontology_apiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('cb_annotation_ontology_api'):
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
                            {'service': 'cb_annotation_ontology_api',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = cb_annotation_ontology_api(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        f = open('/kb/module/test/debug.json')
        testparams = json.load(f)
        self.serviceImpl.add_annotation_ontology_events(self.ctx,testparams)
        
        output = self.serviceImpl.add_annotation_ontology_events(self.ctx,{
            "input_ref" : "Escherichia_coli_K-12_MG1655",
            "input_workspace" : 68056,
            "output_name" : "TestOutput",
            "events" : [{
                #"event_id": "RAST-annotate_genome:SSO:2022-05-11T07:19:52",
                #"original_description": "RAST annotation:2022-05-11T07:19:52",
                "description": "Test",
                "ontology_id": "EC",
                "method": "add_annotation_ontology_events",
                "method_version": "0.0.9",
                "timestamp": "2022-05-11T07:19:52",
                "ontology_terms": {
                    "A0A031WDE4_gene": [
                      {
                        "term": "EC:4.2.1.172",
                        "modelseed_ids": [
                          "MSRXN:rxn42623"
                        ]
                      }
                    ]
                }
            }],
            "output_workspace": 68056,
            "save" : 1
        })
        output = self.serviceImpl.get_annotation_ontology_events(self.ctx,{
            "input_ref" : "TestOutput",
            "input_workspace" : 68056
        })