from urllib.request import urlopen
import time
import json
import os
import re

PS_url = 'https://raw.githubusercontent.com/ModelSEED/PlantSEED/'
PS_tag = 'kbase_release'

# These should be retrieved from the Template data
Template_Compartment_Mapping={'c':'cytosol', 'g':'golgi', 'w':'cellwall',
                              'n':'nucleus', 'r':'endoplasm',
                              'v':'vacuole', 'cv':'vacuole',
                              'd':'plastid', 'cd':'plastid',
                              'm':'mitochondria','cm':'mitochondria',
                              'mj':'mitointer', 'ce':'extracellular',
                              'x':'peroxisome', 'cx':'peroxisome',
                              'e':'extracellular','de':'plastid'}

class FetchPlantSEEDImpl:

    def fetch_reactions(self):

        reactions_data = dict()

        # Load these directly from PlantSEED_Roles.json
        PS_Roles = json.load(urlopen(PS_url+PS_tag+'/Data/PlantSEED_v3/PlantSEED_Roles.json'))

        for entry in PS_Roles:
            if(entry['include'] is False):
                # print(entry['role'])
                continue

            main_class_ss = list()
            main_class = list()
            for metabolic_class in entry['classes']:
                if(len(entry['classes'][metabolic_class].keys())>0):
                    main_class.append(metabolic_class)

                for ss in entry['classes'][metabolic_class].keys():
                    if(ss not in main_class_ss):
                        main_class_ss.append(ss)

            for rxn in entry['reactions']:
                if(rxn not in reactions_data):
                    reactions_data[rxn]={'ecs':[],
                                         'roles':[],
                                         'classes':[],
                                         'subsystems':[],
                                         'compartments':[]}

                if(entry['role'] not in reactions_data[rxn]['roles']):
                    reactions_data[rxn]['roles'].append(entry['role'])

                for mclass in main_class:
                    if(mclass not in reactions_data[rxn]['classes']):
                        reactions_data[rxn]['classes'].append(mclass)

                for subsystem in main_class_ss:
                    if(subsystem not in reactions_data[rxn]['subsystems']):
                        reactions_data[rxn]['subsystems'].append(subsystem)

                for cpt in entry['localization']:
                    if(cpt not in reactions_data[rxn]['compartments']):
                        reactions_data[rxn]['compartments'].append(cpt)

        for rxn in reactions_data:
            for role in reactions_data[rxn]['roles']:
                if('EC' in role):
                    match = re.search(r"\d+\.[\d-]+\.[\d-]+\.[\d-]+", role)
                    if(match is not None):
                        if(match.group(0) not in reactions_data[rxn]['ecs']):
                            reactions_data[rxn]['ecs'].append(match.group(0))

        print("Collected "+str(len(reactions_data))+" core reactions")
        return reactions_data


    def fetch_features(self):

        features_data = dict()

        # Load these directly from PlantSEED_Roles.json
        PS_Roles = json.load(urlopen(PS_url+PS_tag+'/Data/PlantSEED_v3/PlantSEED_Roles.json'))

        for entry in PS_Roles:
        #    if(entry['include'] is False):
        #        continue

            for feature in entry['features']:
                if(feature not in features_data):
                    features_data[feature]= {'roles':[],
                                             'compartments':[]}

                if(entry['role'] not in features_data[feature]['roles']):
                    features_data[feature]['roles'].append(entry['role'])

                for compartment in entry['localization']:
                    for curated_ftr in entry['localization'][compartment]:
                        if(curated_ftr == feature):
                            if(compartment not in features_data[feature]['compartments']):
                                features_data[feature]['compartments'].append(compartment)

        # consolidate and create functions
        # this is the heart of all PlantSEED annotation
        for feature in features_data:

            # compose function as sorted roles
            feature_function = " / ".join( sorted( features_data[feature]['roles'] ) )

            # collect mapped compartments
            feature_compartment_list =  list()
            for compartment in sorted( features_data[feature]['compartments'] ):
                if(compartment not in Template_Compartment_Mapping):
                    print("WARNING, missing compartment: ",compartment)
                mapped_compartment = Template_Compartment_Mapping[compartment]
                feature_compartment_list.append(mapped_compartment)

            # compose sorted list of compartments
            # NB: these are sorted according to single-letter identifier
            feature_compartments = ""
            if(len(feature_compartment_list)>0):
                feature_compartments = " # "+" # ".join( feature_compartment_list )

            full_function = feature_function+feature_compartments
            features_data[feature]['function']=full_function

        return(features_data)

    def log(self,message, prefix_newline=False):
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
        print(('\n' if prefix_newline else '') + time_str + ': ' + message)

    def __init__(self):
        pass

def main():
    plantseed = FetchPlantSEEDImpl()
    plantseed_reactions = plantseed.fetch_reactions()

    plantseed_features = plantseed.fetch_features()

    # Cross-checking plantseed functional annotation with that originally used by OrthoFinder
    PS_Annotation_File = '/Data/PlantSEED_v3/Functional_Annotation/Arabidopsis_Family_Curation.txt'
    PS_Annotation_Path = PS_url+PS_tag+PS_Annotation_File
    with urlopen(PS_Annotation_Path) as fh:
        for line in fh.readlines():
            line = line.decode("utf-8")
            line = line.strip('\r\n')
            
            (family,transcript,curation)=line.split('\t')
            if(curation == "Uncurated"):
                continue
                
            gene = ".".join(transcript.split('.')[0:-1])
            if(curation != plantseed_features[gene]['function']):
                # As of 08/26/21, there is a 100% match so this doesn't trigger
                print(gene,">"+curation+"<",">"+plantseed_features[gene]['function']+"<")

if(__name__ == "__main__"):
    main()
