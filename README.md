# cb_annotation_ontology_api

Provides tooling for adding and retrieving ontology terms from KBase genomes. See the sp

## Current Status

| Branch  | Build                                                              | Coverage                                                                         | LGTM Alerts                                                     |
| ------- | ------------------------------------------------------------------ | -------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| master  | [![KBase SDK Tests](https://github.com/kbaseapps/cb_annotation_ontology_api/workflows/KBase%20SDK%20Tests/badge.svg)](https://github.com/kbaseapps/cb_annotation_ontology_api/actions?query=workflow%3A%22KBase+SDK+Tests%22)  | [![codecov](https://codecov.io/gh/kbaseapps/cb_annotation_ontology_api/branch/master/graph/badge.svg)](https://codecov.io/gh/kbaseapps/cb_annotation_ontology_api)  | [![Total alerts](https://img.shields.io/lgtm/alerts/g/kbaseapps/cb_annotation_ontology_api.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/kbaseapps/cb_annotation_ontology_api/alerts/)  |

## Documentation

This API provides support for loading and retrieving ontology terms to and from objects in the KBase workspace. Currently just the Genome object is supported, but soon AnnotatedMetagenomeAssembly, DNASequenceSet, and ProteinSequenceSet
will all be supported as well. Ontology data-structures can be complex and have many properties they must conform to. This API manages all of these. Also, ontologies can be mapped to ModelSEED reaction IDs to facilitate comparison, and
this API provides that functionality as well. The primary functions for the API are described below, with details outlined in the [module spec file](https://github.com/kbaseapps/cb_annotation_ontology_api/blob/main/cb_annotation_ontology_api.spec). 

## Functions

### get_annotation_ontology_events({
	string input_ref - a reference to the object from which ontology terms should be retrieved
	string input_workspace - name or numerical ID of the workspace containing the object (optional if full ref is provided above)
	list<string> query_events - list of specific events to retrieve (optional - all events are retrieved if empty)
	list<string> query_genes - list of genes to pull data for (optional - all events are retrieved if empty)
	int standardize_modelseed_ids - specify 1 if the ModelSEED IDs should be standardized
}) - this function retrieves the ontology terms from a KBase object

### add_annotation_ontology_events({
	string input_ref - a reference to the object to which ontology terms should be added
	string input_workspace - name or numerical ID of the workspace containing the object (optional if full ref is provided above)
	string output_name - new name to which the modified object should be saved
	string output_workspace - new name or numerical ID of workspace to which object should be saved(optional - input workspace will be used if left blank)
	int clear_existing - set this flag to "1" to clear existing ontology terms from the object (0 is default)
	int overwrite_matching - set this flag to "1" to overwrite existing ontology terms with the same ID in the object (1 is default)
	list<AnnotationOntologyEvent> events - datastructure containing new ontology terms to be written to the object
}) - this function adds new ontology terms to a KBase object


#### AnnotationOntologyEvent structure { 
	string event_id;
	string description;
	string ontology_id;
	string method;
	string method_version;
	string timestamp;
	mapping<string gene_id,string type> feature_types;
	mapping<string gene_id,list<AnnotationOntologyTerm> terms> ontology_terms
}

#### AnnotationOntologyTerm structure {
	string term;
    list<string> modelseed_ids;
    string evidence;
}

## General SDK notes

This is a [KBase](https://kbase.us) module generated by the [KBase Software Development Kit (SDK)](https://github.com/kbase/kb_sdk).

You will need to have the SDK installed to use this module. [Learn more about the SDK and how to use it](https://kbase.github.io/kb_sdk_docs/).

You can also learn more about the apps implemented in this module from its [catalog page](https://narrative.kbase.us/#catalog/modules/cb_annotation_ontology_api) or its [spec file]($module_name.spec).

# Setup and test

Add your KBase developer token to `test_local/test.cfg` and run the following:

```bash
$ make
$ kb-sdk test
```

After making any additional changes to this repo, run `kb-sdk test` again to verify that everything still works.

# Installation from another module

To use this code in another SDK module, call `kb-sdk install cb_annotation_ontology_api` in the other module's root directory.

# Help

You may find the answers to your questions in our [FAQ](https://kbase.github.io/kb_sdk_docs/references/questions_and_answers.html) or [Troubleshooting Guide](https://kbase.github.io/kb_sdk_docs/references/troubleshooting.html).



