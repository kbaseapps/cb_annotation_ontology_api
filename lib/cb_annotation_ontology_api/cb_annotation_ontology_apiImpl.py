# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import json
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService
from cb_annotation_ontology_api.annotation_ontology_api import AnnotationOntologyAPI
#END_HEADER


class cb_annotation_ontology_api:
    '''
    Module Name:
    cb_annotation_ontology_api

    Module Description:
    A KBase module: cb_annotation_ontology_api
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.dfu_client = DataFileUtil(self.config['SDK_CALLBACK_URL'])
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        self.ws_client = workspaceService(config["workspace-url"])
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def get_annotation_ontology_events(self, ctx, params):
        """
        Retrieves annotation ontology events in a standardized form cleaning up inconsistencies in underlying data
        :param params: instance of type "GetAnnotationOntologyEventsParams"
           -> structure: parameter "input_ref" of String, parameter
           "input_workspace" of String, parameter "query_events" of list of
           String, parameter "query_genes" of list of String, parameter
           "standardize_modelseed_ids" of Long
        :returns: instance of type "GetAnnotationOntologyEventsOutput" ->
           structure: parameter "events" of list of type
           "AnnotationOntologyEvent" -> structure: parameter "event_id" of
           String, parameter "description" of String, parameter "ontology_id"
           of String, parameter "method" of String, parameter
           "method_version" of String, parameter "timestamp" of String,
           parameter "feature_types" of mapping from String to String,
           parameter "ontology_terms" of mapping from String to list of type
           "AnnotationOntologyTerm" -> structure: parameter "term" of String,
           parameter "modelseed_ids" of list of String, parameter "evidence"
           of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN get_annotation_ontology_events
        anno_api = AnnotationOntologyAPI(self.config,self.ws_client,self.dfu_client)
        output = anno_api.get_annotation_ontology_events(params)
        with open(self.shared_folder+'/debug.json', 'w') as outfile:
            json.dump(output, outfile,indent=4)
        #END get_annotation_ontology_events

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method get_annotation_ontology_events return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def add_annotation_ontology_events(self, ctx, params):
        """
        Adds a new annotation ontology event to a genome or AMA
        :param params: instance of type "AddAnnotationOntologyEventsParams"
           -> structure: parameter "input_ref" of String, parameter
           "input_workspace" of String, parameter "output_name" of String,
           parameter "output_workspace" of String, parameter "clear_existing"
           of Long, parameter "overwrite_matching" of Long, parameter
           "events" of list of type "AnnotationOntologyEvent" -> structure:
           parameter "event_id" of String, parameter "description" of String,
           parameter "ontology_id" of String, parameter "method" of String,
           parameter "method_version" of String, parameter "timestamp" of
           String, parameter "feature_types" of mapping from String to
           String, parameter "ontology_terms" of mapping from String to list
           of type "AnnotationOntologyTerm" -> structure: parameter "term" of
           String, parameter "modelseed_ids" of list of String, parameter
           "evidence" of String
        :returns: instance of type "AddAnnotationOntologyEventsOutput" ->
           structure: parameter "output_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN add_annotation_ontology_events
        anno_api = AnnotationOntologyAPI(self.config,self.ws_client,self.dfu_client)
        output = anno_api.add_annotation_ontology_events(params)
        #END add_annotation_ontology_events

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method add_annotation_ontology_events return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
