/*
A KBase module: cb_annotation_ontology_api
*/

module cb_annotation_ontology_api {
    typedef structure {
    	string term;
    	list<string> modelseed_ids;
    	string evidence;
    } AnnotationOntologyTerm;
    
    typedef structure {
		string event_id;
		string description;
		string ontology_id;
		string method;
		string method_version;
		string timestamp;
		mapping<string gene_id,string type> feature_types;
		mapping<string gene_id,list<AnnotationOntologyTerm> terms> ontology_terms;
	} AnnotationOntologyEvent;
    
    typedef structure {
		string input_ref;
		string input_workspace;
		list<string> query_events;
		list<string> query_genes;
		int standardize_modelseed_ids;
    } GetAnnotationOntologyEventsParams;
    
    typedef structure {
		list<AnnotationOntologyEvent> events;
    } GetAnnotationOntologyEventsOutput;
    
    /*
        Retrieves annotation ontology events in a standardized form cleaning up inconsistencies in underlying data
    */
    funcdef get_annotation_ontology_events(GetAnnotationOntologyEventsParams params) returns (GetAnnotationOntologyEventsOutput output) authentication optional;
	
	typedef structure {
		string input_ref;
		string input_workspace;
		string output_name;
		string output_workspace;
		int clear_existing;
		int overwrite_matching;
		list<AnnotationOntologyEvent> events;
    } AddAnnotationOntologyEventsParams;
    
    typedef structure {
		string output_ref;
    } AddAnnotationOntologyEventsOutput;
    
    /*
        Adds a new annotation ontology event to a genome or AMA
    */
	funcdef add_annotation_ontology_events(AddAnnotationOntologyEventsParams params) returns (AddAnnotationOntologyEventsOutput output) authentication optional;

};
