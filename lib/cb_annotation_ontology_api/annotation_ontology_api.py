import json
import os
import re

# silence whining
import requests
import hashlib
requests.packages.urllib3.disable_warnings()

source_hash = {
    "MetaCyc" : "META",
    "KEGG" : "RO",
    "BiGG" : "BIGG",
    "rhea" : "RHEA"
}

ontology_translation = {
    "KEGGKO" : "KO",
    "KEGGRO" : "RO",
    "METACYC" : "META",
    "SEED" : "SSO",
    "TCDB" : "TC",
    "MODELSEED" : "MSRXN"
}

ontology_hash = {
    "KO" : 1,
    "EC" : 1,
    "SSO" : 1,
    "RO" : 1,
    "META" : 1,
    "MSRXN" : 1,
    "MSCPD" : 1,
    "MSCPX" : 1,
    "BIGG" : 1,
    "BIGGCPD" : 1,
    "GO" : 1,
    "TC" : 1,
    "RHEA" : 1
};

class AnnotationOntologyAPI:
    def __init__(self,config,ws_client = None, dfu_client = None,gfu_client = None):
        self.ws_client = ws_client
        self.dfu_client = dfu_client
        self.gfu_client = gfu_client
        self.alias_hash = {}
        self.term_names = {}
        self.config = config
    
    def process_workspace_identifiers(self,id_or_ref, workspace=None):
        """
        IDs should always be processed through this function so we can interchangeably use
        refs, IDs, and names for workspaces and objects
        """
        objspec = {}
        if workspace is None or len(id_or_ref.split("/")) > 1:
            objspec["ref"] = id_or_ref
        else:
            if isinstance(workspace, int):
                objspec['wsid'] = workspace
            else:
                objspec['workspace'] = workspace
            if isinstance(id_or_ref, int):
                objspec['objid'] = id_or_ref
            else:
                objspec['name'] = id_or_ref
                
        #print("Object spec:")
        #for key in objspec:
        #    print(key+"\t"+objspec[key])
        return objspec
    
    def get_alias_hash(self,namespace):
        if "MSRXN" not in self.alias_hash:
            filename = self.config["data_directory"]+"/msrxn_hash.json"
            with open(filename) as json_file:
                self.alias_hash["MSRXN"] = json.load(json_file)
        if namespace not in self.alias_hash:
            self.alias_hash[namespace] = {}
            if namespace == "EC":
                filename = self.config["data_directory"]+"/EC_translation.tsv"
                data = ""
                with open(filename, 'r') as file:
                    data = file.read()
                lines = data.split("\n")
                lines.pop(0)
                for line in lines:
                    items = line.split("\t")
                    if len(items) >= 2:
                        modelseed = "MSRXN:"+items[0]
                        if modelseed in self.alias_hash["MSRXN"]:
                            modelseed = self.alias_hash["MSRXN"][modelseed][0]
                        if items[1] not in self.alias_hash["EC"]:
                            self.alias_hash["EC"]["EC:"+items[1]] = []
                        self.alias_hash["EC"]["EC:"+items[1]].append(modelseed)
            elif namespace == "META" or namespace == "RO" or namespace == "BIGG" or namespace == "RHEA":
                filename = self.config["data_directory"]+"/ModelSEED_Reaction_Aliases.txt"
                data = ""
                with open(filename, 'r') as file:
                    data = file.read()
                lines = data.split("\n")
                for line in lines:
                    items = line.split("\t")
                    if len(items) >= 3:
                        modelseed = "MSRXN:"+items[0]
                        if modelseed in self.alias_hash["MSRXN"]:
                            modelseed = self.alias_hash["MSRXN"][modelseed][0]
                        source = None
                        if items[2] in source_hash:
                            source = source_hash[items[2]]
                        if source != None:
                            if source not in self.alias_hash:
                                self.alias_hash[source] = {}
                            if items[1] not in self.alias_hash[source]:
                                self.alias_hash[source][source+":"+items[1]] = []
                            self.alias_hash[source][source+":"+items[1]].append(modelseed)
            elif namespace == "KO":
                filename = self.config["data_directory"]+"/kegg_95_0_ko_seed.tsv"
                data = ""
                with open(filename, 'r') as file:
                    data = file.read()
                lines = data.split("\n")
                lines.pop(0)
                for line in lines:
                    items = line.split("\t")
                    if len(items) >= 2:
                        if items[0] not in self.alias_hash["KO"]:
                            self.alias_hash["KO"]["KO:"+items[0]] = []
                        modelseed_ids = items[1].split(";")
                        id_hash = {}
                        for modelseed in modelseed_ids:
                            modelseed = "MSRXN:"+modelseed
                            if modelseed in self.alias_hash["MSRXN"]:
                                modelseed = self.alias_hash["MSRXN"][modelseed][0]
                            if modelseed not in id_hash:
                                self.alias_hash["KO"]["KO:"+items[0]].append(modelseed)
                            id_hash[modelseed] = 1
            elif namespace == "SSO":
                sso_template = dict()
                filename = self.config["data_directory"]+"/SSO_reactions.json"
                with open(filename) as json_file:
                    sso_template = json.load(json_file)
                for sso in sso_template:
                    id_hash = {}
                    for modelseed in sso_template[sso]:
                        modelseed = "MSRXN:"+modelseed
                        if modelseed in self.alias_hash["MSRXN"]:
                            modelseed = self.alias_hash["MSRXN"][modelseed][0]
                        if modelseed not in id_hash:
                            if sso not in self.alias_hash["SSO"]:
                                self.alias_hash["SSO"][sso] = []
                            self.alias_hash["SSO"][sso].append(modelseed)
                        id_hash[modelseed] = 1             
            elif namespace == "GO":
                go_translation = dict()
                filename = self.config["data_directory"]+"/GO_ontology_translation.json"
                with open(filename) as json_file:
                    go_translation = json.load(json_file)
                for term in go_translation["translation"]:
                    adjusted_term = "GO:"+term
                    if "equiv_terms" in go_translation["translation"][term]:
                        id_hash = {}
                        for rxn_data in go_translation["translation"][term]["equiv_terms"]:
                            if rxn_data["equiv_term"] != None:
                                modelseed = "MSRXN:"+rxn_data["equiv_term"]
                                if modelseed in self.alias_hash["MSRXN"]:
                                    modelseed = self.alias_hash["MSRXN"][modelseed][0]
                                if adjusted_term not in self.alias_hash["GO"]:
                                    self.alias_hash["GO"][adjusted_term] = []
                                if modelseed not in id_hash:
                                    self.alias_hash["GO"][adjusted_term].append(modelseed)
                                id_hash[modelseed] = 1            
        return self.alias_hash[namespace]
                
    def translate_term_to_modelseed(self,term):
        namespace = term.split(":").pop(0)
        if namespace == "MSRXN":
            if term not in self.get_alias_hash(namespace):
                return [term]
            else:
                return self.get_alias_hash(namespace)[term]
        elif term not in self.get_alias_hash(namespace):
            return []
        else:
            return self.get_alias_hash(namespace)[term]
        
    def get_annotation_ontology_events(self,params):
        #Building query hash
        event_query = None
        if "query_events" in params and not params["query_events"] == None:
            for event in params["query_events"]:
                event_query[event] = 1
        gene_query = None
        if "query_genes" in params and not params["query_genes"] == None:
            for gene in params["query_genes"]:
                gene_query[gene] = 1
        #Pull the object from the workspace is necessary
        if "object" not in params:
            res = None
            if "input_workspace" not in params:
                res = self.ws_client.get_objects2({"objects": [self.process_workspace_identifiers(params["input_ref"], None)]})
            else: 
                res = self.ws_client.get_objects2({"objects": [self.process_workspace_identifiers(params["input_ref"], params["input_workspace"])]})
            params["object"] = res["data"][0]["data"]
            params["type"] = res["data"][0]["info"][2]
        #Get the feature data
        features = []
        types = {}
        if "features" in params["object"]:
            features.extend(params["object"]["features"])
            for ftr in params["object"]["features"]:
                types[ftr["id"]] = "gene"
        if "cdss" in params["object"]:
            features.extend(params["object"]["cdss"])
            for ftr in params["object"]["cdss"]:
                types[ftr["id"]] = "cds"
        if "mrnas" in params["object"]:
            features.extend(params["object"]["mrnas"])
            for ftr in params["object"]["mrnas"]:
                types[ftr["id"]] = "mrna"
        if "non_coding_features" in params["object"]:
            features.extend(params["object"]["non_coding_features"])
            for ftr in params["object"]["non_coding_features"]:
                types[ftr["id"]] = "noncoding"
        elif "features_handle_ref" in params["object"]:
            shock_output = self.dfu_client.shock_to_file({
                "handle_id" : params["object"]["features_handle_ref"],
                "file_path" : self.config["scratch"]
            })
            os.system("gunzip --force ".shock_output["file_path"])
            shock_output["file_path"] = shock_output["file_path"][0:-3]
            with open(shock_output["file_path"]) as json_file:
                features = json.load(json_file)
            for ftr in features:
                types[ftr["id"]] = "gene"
        output = {"events" : [],"feature_types" : {}}
        if "ontology_events" in params["object"]:
            events_array = []
            for event in params["object"]["ontology_events"]:
                if "event_id" not in event:
                    event["event_id"] = event["method"]+":"+event["method_version"]+":"+event["id"]+":"+event["timestamp"]
                old_description = None
                if "description" in event:
                    old_description = event["description"]
                    if event["description"][-1*len(event["timestamp"]):] != event["timestamp"]:
                        event["description"] = event["description"]+":"+event["timestamp"]
                else:
                    event["description"] = event["method"]+":"+event["method_version"]+":"+event["id"]+":"+event["timestamp"]
                newevent = {
                    "event_id" : event["event_id"],
                    "original_description" : old_description,
                    "description" : event["description"],
                    "ontology_id" : event["id"].upper(),
                    "method" : event["method"],
                    "method_version" : event["method_version"],
                    "timestamp" : event["timestamp"],
                    "ontology_terms" : {}
                }
                if newevent["ontology_id"] not in ontology_hash and newevent["ontology_id"] in ontology_translation:
                    newevent["ontology_id"] = ontology_translation[newevent["ontology_id"]]
                events_array.append(newevent)
                if event_query == None or id in event_query:
                    output["events"].append(newevent)
            for feature in features:
                if gene_query == None or feature["id"] in gene_query:
                    if "ontology_terms" in feature:
                        for tag in feature["ontology_terms"]:
                            original_tag = tag
                            tag = tag.upper()
                            if tag not in ontology_hash and tag in ontology_translation:
                                tag = ontology_translation[tag]
                            if tag in ontology_hash:
                                for term in feature["ontology_terms"][original_tag]:
                                    original_term = term
                                    array = term.split(":")
                                    if len(array) == 1:
                                        term = tag+":"+array[0]
                                    else:
                                        if array[0].upper() == original_tag.upper() or array[0].upper() == tag:
                                            array[0] = tag
                                            term = ":".join(array)
                                        else:
                                            term = tag+":"+":".join(array)
                                    modelseed_ids = self.translate_term_to_modelseed(term)
                                    termhash = {}
                                    for event_index in feature["ontology_terms"][original_tag][original_term]:
                                        if feature["id"] not in events_array[event_index]["ontology_terms"]:
                                            output["feature_types"][feature["id"]] = types[feature["id"]]
                                            events_array[event_index]["ontology_terms"][feature["id"]] = []
                                        if term not in termhash:
                                            termhash[term] = {}
                                        termhash[term][event_index] = 1
                                    for term in termhash:
                                        for event_index in termhash[term]:
                                            termdata = {"term" : term}
                                            if len(modelseed_ids) > 0:
                                                termdata["modelseed_ids"] = modelseed_ids
                                            if "ontology_evidence" in feature:
                                                if original_term in feature["ontology_evidence"]:
                                                    if event_index in feature["ontology_evidence"][original_term]:
                                                        termdata["evidence"] = feature["ontology_evidence"][original_term][event_index]
                                            events_array[event_index]["ontology_terms"][feature["id"]].append(termdata)
        return output
    
    def add_annotation_ontology_events(self,params):
        #Pull the object from the workspace is necessary
        ref = None
        if "object" not in params or params["object"] == None:
            if "input_workspace" not in params:
                res = self.ws_client.get_objects2({"objects": [self.process_workspace_identifiers(params["input_ref"], None)]})
            else: 
                res = self.ws_client.get_objects2({"objects": [self.process_workspace_identifiers(params["input_ref"], params["input_workspace"])]})
            params["object"] = res["data"][0]["data"]
            params["type"] = res["data"][0]["info"][2]
            ref = str(res["data"][0]["info"][6])+"/"+str(res["data"][0]["info"][0])+"/"+str(res["data"][0]["info"][4])
        output = {
            "ftrs_not_found" : [],"ftrs_found" : 0,"terms_not_found" : []
        }
        #Pulling existing ontology so we can standardize and check for matches
        ontologies_present = {}
        events = self.get_annotation_ontology_events(params)["events"]
        if "clear_existing" in params and params["clear_existing"] == 1: 
            events = []
        #Scrolling through new events, stadardizing, and checking for matches
        new_events = {}
        for event in params["events"]:
            if "ontology_id" not in event:
                event["ontology_id"] = event["id"]
            event["ontology_id"] = event["ontology_id"].upper()
            if event["ontology_id"] in ontology_translation:
                event["ontology_id"] = ontology_translation[event["ontology_id"]]
            event["id"] = event["ontology_id"]
            new_events[event["id"]] = 1
            #Creating description
            if "event_id" not in event:
                event["event_id"] = event["method"]+":"+event["ontology_id"]+":"+event["timestamp"]
            if "description" not in event:
                event["description"] = event["method"]+":"+event["method_version"]+":"+event["ontology_id"]+":"+event["timestamp"]
            elif event["description"][-1*len(event["timestamp"]):] != event["timestamp"]:
                event["description"] = event["description"]+":"+event["timestamp"]
            index = 0
            match = 0
            for existing_event in events:
                #If an existing event has a matching event ID, we overwrite it
                if existing_event["event_id"] == event["event_id"]:
                    match = 1
                    if "overwrite_matching" in params and params["overwrite_matching"] == 1:
                        events[index] = event
                index += 1
            if match == 0:
                events.append(event)
        #Filling feature hash with all feature types which should all have unique ids
        alias_hash = {}
        feature_hash = {}
        lc_feature_hash = {}
        lc_alias_hash = {}
        terms_not_found = {}
        feature_found_hash = {}
        feature_types = ["features","cdss","mrnas","non_coding_features"]
        for currtype in feature_types:
            if currtype in params["object"]:
                to_remove = []
                for ftr in params["object"][currtype]: 
                    ftr["ontology_terms"] = {}
                    if currtype == "features" and "protein_translation" not in ftr:
                        if "non_coding_features" not in params["object"]:
                            params["object"]["non_coding_features"] = []
                        params["object"]["non_coding_features"].append(ftr)
                        to_remove.append(ftr)
                    else:
                        self.upgrade_feature(ftr,currtype)
                        feature_hash[ftr["id"]] = ftr
                        lc_feature_hash[ftr["id"].lower()] = ftr
                        self.process_feature_aliases(ftr,alias_hash,lc_alias_hash)
                for item in to_remove:
                    params["object"][currtype].remove(item)
        if "features_handle_ref" in params["object"]:
            if "feature_object" not in params:
                shock_output = self.dfu_client.shock_to_file({
                    "handle_id" : params["object"]["features_handle_ref"],
                    "file_path" : self.config["scratch"]
                })
                os.system("gunzip --force ".shock_output["file_path"])
                shock_output["file_path"] = shock_output["file_path"][0:-3]
                with open(shock_output["file_path"]) as json_file:
                    params["feature_object"] = json.load(json_file)
        if "feature_object" in params:
            for ftr in params["feature_object"]:
                feature_hash[ftr["id"]] = ftr
                lc_feature_hash[ftr["id"].lower()] = ftr
                self.process_feature_aliases(ftr,alias_hash,lc_alias_hash)
        #Adding events
        params["object"]["ontology_events"] = []
        for event in events:
            new_event = {
                "description" : event["description"],
                "id" : event["ontology_id"],
                "event_id" : event["event_id"],
                "ontology_id" : event["ontology_id"],
                "method" : event["method"],
                "method_version" : event["method_version"],
                "timestamp" : event["timestamp"]
            }
            if "ontology_events" not in params["object"]:
                params["object"]["ontology_events"] = []
            event_index = len(params["object"]["ontology_events"])
            params["object"]["ontology_events"].append(new_event)
            for currgene in event["ontology_terms"]:
                genes = []
                if currgene in feature_hash:
                    if new_event["id"] in new_events:
                        feature_found_hash[currgene] = 1
                    genes = [currgene]
                elif currgene in alias_hash:
                    if new_event["id"] in new_events:
                        feature_found_hash[currgene] = 1
                    genes = alias_hash[currgene]
                elif currgene.lower() in lc_feature_hash:
                    if new_event["id"] in new_events:
                        feature_found_hash[currgene] = 1
                    genes = [lc_feature_hash[currgene.lower()]["id"]]
                elif currgene.lower() in lc_alias_hash:
                    if new_event["id"] in new_events:
                        feature_found_hash[currgene] = 1
                    genes = lc_alias_hash[currgene.lower()]
                else:
                    output["ftrs_not_found"].append(currgene)
                for gene in genes:
                    if gene in feature_hash:
                        feature = feature_hash[gene]
                        if "ontology_terms" not in feature:
                            feature["ontology_terms"] = {}
                        if new_event["id"] not in feature["ontology_terms"]:
                            feature["ontology_terms"][new_event["id"]] = {}
                        for term in event["ontology_terms"][currgene]:
                            if term["term"].split(":")[0] != new_event["id"]:
                                term["term"] = new_event["id"]+":"+term["term"]
                            #If this is a SEED role, translate to an SSO
                            if new_event["id"] == "SSO" and re.search('^SSO:\d+$', term["term"]) == None:
                                term["term"] = re.sub("^SSO:","",term["term"])
                                terms = re.split("\s*;\s+|\s+[\@\/]\s+",term["term"])
                                first = 1
                                for subterm in terms:
                                    if first == 1:
                                        #Only the first term completes the rest of this code
                                        term["term"] = self.translate_rast_function_to_sso(subterm)
                                        first = 0
                                    else:
                                        #Sub sterms need to be added independently
                                        subterm = self.translate_rast_function_to_sso(subterm)
                                        if subterm != None:
                                            if subterm not in feature["ontology_terms"][new_event["id"]]:
                                                feature["ontology_terms"][new_event["id"]][subterm] = []
                                            feature["ontology_terms"][new_event["id"]][subterm].append(event_index)
                                            if new_event["id"] not in ontologies_present:
                                                ontologies_present[new_event["id"]] = {}
                                            ontologies_present[new_event["id"]][subterm] = self.get_term_name(new_event["id"],subterm)                        
                                            if ontologies_present[new_event["id"]][subterm] == "Unknown":
                                                terms_not_found[subterm] = 1
                            if term["term"] == None:
                                continue
                            if term["term"] not in feature["ontology_terms"][new_event["id"]]:
                                feature["ontology_terms"][new_event["id"]][term["term"]] = []
                            feature["ontology_terms"][new_event["id"]][term["term"]].append(event_index)
                            if new_event["id"] not in ontologies_present:
                                ontologies_present[new_event["id"]] = {}
                            ontologies_present[new_event["id"]][term["term"]] = self.get_term_name(new_event["id"],term["term"])
                            if ontologies_present[new_event["id"]][term["term"]] == "Unknown":
                                terms_not_found[term["term"]] = 1
                            if "evidence" in term:
                                if "ontology_evidence" not in feature:
                                    feature["ontology_evidence"] = {}
                                if term["term"] not in feature["ontology_evidence"]:
                                    feature["ontology_evidence"][term["term"]] = {}
                                feature["ontology_evidence"][term["term"]][event_index] = term["evidence"]
        output["ftrs_found"] = len(feature_found_hash)
        for term in terms_not_found:
            output["terms_not_found"].append(term)
        params["object"]["ontologies_present"] = ontologies_present
        #Saving object if requested but not if it's an AMA
        if params["save"] == 1:
            #Setting provenance
            provenance_params = {}
            for key in params:
                if not key == "object" and not key == "events" and not "feature_object":
                    provenance_params[key] = params[key]            
            provenance = [{
                'description': 'A function that adds ontology terms to a genome or metagenome',
                'input_ws_objects': [],
                'method': 'add_annotation_ontology_events',
                'method_params': [provenance_params],
                'service': 'annotation_ontology_api',
                'service_ver': 1,
            }]
            #If a metagenome, saving features
            params["type"] = "KBaseGenomes.Genome"
            if "feature_object" in params:
                params["type"] = "KBaseMetagenomes.AnnotatedMetagenomeAssembly"
                json_file_path = self.config["scratch"]+params["object"]["name"]+"_features.json"
                with open(json_file_path, 'w') as fid:
                    json.dump(params["feature_object"], fid)
                json_to_shock = self.dfu_client.file_to_shock(
                    {'file_path': json_file_path, 'make_handle': 1, 'pack': 'gzip'}
                )
                # Resetting feature file handle o new value
                params["object"]['features_handle_ref'] = json_to_shock['handle']['hid']
                # Remove json file to avoid disk overload
                os.remove(json_file_path)
            # Removing genbank handle ref because this breaks saving
            params["object"].pop('genbank_handle_ref', None)
            #Adding missing fields in genome
            if params["type"] == "KBaseGenomes.Genome":
                self.check_genome(params["object"],ref)
            # Saving genome/metagenome object to workspace
            gfu_param = {
                "name" : params["output_name"],
                "data" : params["object"],
                "upgrade" : 1,
                "provenance" : params["provenance"],
                "hidden" : 0
            }
            if isinstance(params["output_workspace"], int):
                gfu_param['id'] = params["output_workspace"]
            else:
                gfu_param['workspace'] = params["output_workspace"]
            save_output = self.gfu_client.save_one_genome(gfu_param);
            output["output_ref"] = str(save_output["info"][6])+"/"+str(save_output["info"][0])+"/"+str(save_output["info"][4])
            output["output_name"] = str(save_output["info"][1])
        else:            
            #Returning object if save not requested
            output["object"] = params["object"]
            output["type"] = params["type"]
            if "feature_object" in params:
                output["feature_object"] = params["feature_object"]
        return output
    
    def process_feature_aliases(self,ftr,alias_hash,lc_alias_hash):
        if "aliases" in ftr:
            for alias in ftr["aliases"]:    
                if not isinstance(alias, str):
                    alias = alias[1]
                if alias not in alias_hash:
                    alias_hash[alias] = []
                alias_hash[alias].append(ftr["id"])
                if alias.lower() not in lc_alias_hash:
                    lc_alias_hash[alias.lower()] = []
                lc_alias_hash[alias.lower()].append(ftr["id"])
        if "db_xrefs" in ftr:
            for alias in ftr["db_xrefs"]:
                if alias[1] not in alias_hash:
                    alias_hash[alias[1]] = []
                alias_hash[alias[1]].append(ftr["id"])
                if alias[1].lower() not in lc_alias_hash:
                    lc_alias_hash[alias[1].lower()] = []
                lc_alias_hash[alias[1].lower()].append(ftr["id"])
    
    def upgrade_feature(self,ftr,type):
        if "function" in ftr:
            ftr["functions"] = re.split("\s*;\s+|\s+[\@\/]\s+",ftr["function"])
            del ftr["function"]
            #Clearing old ontology terms rather than attempting to translate them
            ftr["ontology_terms"] = {}
        if type == "features" and "cdss" not in ftr:
            ftr["cdss"] = []
        if "dna_sequence_length" not in ftr:
            ftr["dna_sequence_length"] = ftr["location"][0][3]
        if "aliases" in ftr:
            for i in range(0, len(ftr["aliases"])):
                if isinstance(ftr["aliases"][i], str):
                    array = ftr["aliases"][i].split(":")
                    if len(array) > 1:
                        ftr["aliases"][i] = array
                    else:
                        ftr["aliases"][i] = ["Unknown",ftr["aliases"][i]]
        if type == "cdss" and "protein_md5" not in ftr:
            ftr["protein_md5"] = hashlib.md5(ftr["protein_translation"].encode()).hexdigest()
        if "md5" not in ftr:
            if "dna_sequence" in ftr:
                ftr["md5"] = hashlib.md5(ftr["dna_sequence"].encode()).hexdigest()
            elif "protein_translation" in ftr:
                ftr["md5"] = hashlib.md5(ftr["protein_translation"].encode()).hexdigest()
            else:
                ftr["md5"] = ""
        
    def check_genome(self,genome,ref = None):
        if "gc_content" in genome and isinstance(genome["gc_content"],str):
            genome["gc_content"] = float(genome["gc_content"])
        if "md5" not in genome:
            genome["md5"] = ""
        if "molecule_type" not in genome:
            genome["molecule_type"] = "DNA"
        if "genome_tiers" not in genome:
            genome["genome_tiers"] = ["ExternalDB","User"]
        if "feature_counts" not in genome:
            genome["feature_counts"] = {"CDS":0,"gene":0,"non-protein_encoding_gene":0}
            if "cdss" in genome:
                genome["feature_counts"]["CDS"] = len(genome["cdss"])
            if "features" in genome:
                genome["feature_counts"]["gene"] = len(genome["features"])
            if "non_coding_features" in genome:
                genome["feature_counts"]["non-protein_encoding_gene"] = len(genome["non_coding_features"])
        if "cdss" not in genome:
            genome["cdss"] = []
        if "assembly_ref" in genome and not ref == None:
            genome["assembly_ref"] = ref+";"+genome["assembly_ref"]
        
    def convert_role_to_searchrole(self,term):
        term = term.lower()
        term = re.sub("\s","",term)
        term = re.sub("[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]*","",term)
        term = re.sub("\#.*$","",term)
        term = re.sub("\(ec:*\)","",term)
        term = re.sub("[\(\)\[\],-]","",term)
        return term
    
    def translate_rast_function_to_sso(self,term):
        #Stripping out SSO prefix if it's present
        term = re.sub("^SSO:","",term)
        term = self.convert_role_to_searchrole(term)
        #Checking for SSO translation file
        if "SEED_ROLE" not in self.alias_hash:
            self.alias_hash["SEED_ROLE"] = {}
            sso_ontology = dict()
            with open(self.config["data_directory"]+"/SSO_dictionary.json") as json_file:
                sso_ontology = json.load(json_file)
            for term in sso_ontology["term_hash"]:
                name = self.convert_role_to_searchrole(sso_ontology["term_hash"][term]["name"])
                self.alias_hash["SEED_ROLE"][name] = term
            
        #Translating
        if term in self.alias_hash["SEED_ROLE"]:
            return self.alias_hash["SEED_ROLE"][term]
        else:
            return None
    
    def get_term_name(self,type,term):
        if type not in self.term_names:
            self.term_names[type] = {}
            if type == "SSO" or type == "EC" or type == "TC" or type == "META" or type == "RO" or type == "KO" or type == "GO":
                with open(self.config["data_directory"]+"/"+type+"_dictionary.json") as json_file:
                    ontology = json.load(json_file)
                    for term in ontology["term_hash"]:
                        self.term_names[type][term] = ontology["term_hash"][term]["name"]
        if term not in self.term_names[type]:
            return "Unknown"
        return self.term_names[type][term]
        