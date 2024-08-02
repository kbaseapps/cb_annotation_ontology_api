# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import json
import sys
sys.path.append("/deps/KBBaseModules/")
from cb_annotation_ontology_api.annotation_ontology_api import AnnotationOntologyModule
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
    VERSION = "1.0.0"
    GIT_URL = "https://github.com/kbaseapps/cb_annotation_ontology_api.git"
    GIT_COMMIT_HASH = "26202f93b10592c75a5aaca79d447a39f0af4923"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = os.environ['KB_AUTH_TOKEN']
        config["version"] = self.VERSION
        self.anno_api = AnnotationOntologyModule("cb_annotation_ontology_api",config,"/kb/module",None,self.token,callback=self.callback_url)
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
           parameter "modelseed_ids" of list of String, parameter
           "evidence_only" of Long, parameter "evidence" of type "Evidence"
           -> structure: parameter "reference" of tuple of size 2: parameter
           "entity_type" of String, parameter "ref_entity" of String,
           parameter "scores" of mapping from String to Double
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN get_annotation_ontology_events
        output = self.anno_api.get_annotation_ontology_events(params)
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
           "evidence_only" of Long, parameter "evidence" of type "Evidence"
           -> structure: parameter "reference" of tuple of size 2: parameter
           "entity_type" of String, parameter "ref_entity" of String,
           parameter "scores" of mapping from String to Double
        :returns: instance of type "AddAnnotationOntologyEventsOutput" ->
           structure: parameter "output_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN add_annotation_ontology_events
        output = self.anno_api.add_annotation_ontology_events(params)
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
