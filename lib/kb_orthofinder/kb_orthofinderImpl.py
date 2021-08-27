# -*- coding: utf-8 -*-
#BEGIN_HEADER

import glob
import os
import re
import uuid
import shutil
import subprocess
import itertools
import time
from urllib.request import urlopen

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService

from kb_orthofinder.core.generate_table_impl import GenerateTableImpl
from kb_orthofinder.core.generate_figure_impl import GenerateFigureImpl
from kb_orthofinder.GenerateFigure import GenerateFigure

from bokeh.plotting import output_file, save

from kb_orthofinder.core.fetch_plantseed_impl import FetchPlantSEEDImpl

#END_HEADER

class kb_orthofinder:
    '''
    Module Name:
    kb_orthofinder

    Module Description:
    A KBase module: kb_orthofinder
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.2"
    GIT_URL = "git@github.com:kbaseapps/kb_orthofinder.git"
    GIT_COMMIT_HASH = "4210203500471209c25c62b544dc0077da862142"

    #BEGIN_CLASS_HEADER

    def generate_figure(self, params):

        Data = dict()
        Reference_Results_File = "Reference_Phytozome_Threshold.txt"
        with open(os.path.join("/kb/module/data",Reference_Results_File)) as reference_results_handle:
            for line in reference_results_handle.readlines():
                line=line.strip()
                (x,group,y)=line.split('\t')
                y=float(y)
                x=float(x)

                if(group not in Data):
                    Data[group]=list()

                Data[group].append((x,y))
        fig_gen = GenerateFigure(Data)

        uuid_string = str(uuid.uuid4())
        figure_data_file_path=os.path.join(self.scratch,uuid_string)
        os.mkdir(figure_data_file_path)
        fig_gen.generate_figure(figure_path=figure_data_file_path,data_point=params)
        return figure_data_file_path

    def log(self,message, prefix_newline=False):
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
        print(('\n' if prefix_newline else '') + time_str + ': ' + message)

    def compute_clusters(self, cluster):
        features = sorted(cluster.keys())
        clustered_sequences=dict()
        for i in range(len(features)-1):
            for j in range(1,len(features)):
                if(i>=j):
                    continue

                Spp1=features[i].split("||")[0]
                Spp2=features[j].split("||")[0]
                if(Spp1 == features[i]):
                    Spp1 = "_".join(features[i].split("_")[0:-1])
                if(Spp2 == features[j]):
                    Spp2 = "_".join(features[j].split("_")[0:-1])

                if(Spp1 == Spp2):
                    continue

                Seq1 = cluster[features[i]]
                Seq2 = cluster[features[j]]

                AA_Match=0.0;
                for k in range(len(Seq1)):
                    if(Seq1[k] == "-" or Seq2[k] == "-"):
                        continue

                    if(Seq1[k] == Seq2[k]):
                        AA_Match+=1.0

                Seq1=Seq1.replace('-','')
                Seq2=Seq2.replace('-','')

                ID1=AA_Match/len(Seq1)
                ID2=AA_Match/len(Seq2)
                Avg = "{0:.2f}".format((ID1+ID2)/2.0)

                if("Athaliana" in Spp1):
                    if(features[j] not in clustered_sequences):
                        clustered_sequences[features[j]]=dict()
                    clustered_sequences[features[j]][features[i]]=Avg
                elif("Athaliana" in Spp2):
                    if(features[i] not in clustered_sequences):
                        clustered_sequences[features[i]]=dict()
                    clustered_sequences[features[i]][features[j]]=Avg

        return clustered_sequences

    def propagate_annotation(self, cluster, threshold, plantseed_curation):

        top_orthologs = dict()
        for ortholog in (cluster.keys()):
            
            if(cluster[ortholog] not in top_orthologs):
                top_orthologs[cluster[ortholog]]=dict()
            top_orthologs[cluster[ortholog]][ortholog]=1

        top_ortholog_seqid="0.00"
        for seq_id in (top_orthologs.keys()):
            if(float(seq_id)>float(top_ortholog_seqid)):
                top_ortholog_seqid=seq_id

        #Key to controlling the propagation of annotation
        #Is the sequence identity of the highest-scoring ortholog good enough?
        if(float(top_ortholog_seqid) < threshold):
            return (None,"Unannotated",0.0)

        top_ortholog=""
        if(len(top_orthologs[top_ortholog_seqid])==1):
            top_ortholog = list(top_orthologs[top_ortholog_seqid].keys())[0]
        else:
            #There are multiple Arabidopsis orthologs that have the same level of seq. id
            #Collect the actual functions and see if there's multiple functions

            Multi_Functions = dict()
            for ortholog in (top_orthologs[top_ortholog_seqid].keys()):
                function="Unannotated"
                ortholog_gene = ".".join(ortholog.split('.')[0:-1])
                if(ortholog_gene not in plantseed_curation):
                    plantseed_curation[ortholog_gene]={'function':function}
                else:
                    function = plantseed_curation[ortholog_gene]['function']
                Multi_Functions[function]=1

            #Rules are:
            # 0) if 1 function, arbitrarily pick ath ortholog
            # 1) if 2 functions, and one is "Unannotated", 
            #    prioritize annotated function
            # 2) if 2 functions, and none is "Unannotated" or > 2 functions, 
            #    its ambiguous, return without doing anything

            if(len(Multi_Functions.keys())==1):
                top_ortholog = list(top_orthologs[top_ortholog_seqid].keys())[0]
            elif(len(Multi_Functions.keys())==2 and "Unannotated" in Multi_Functions.keys()):
                for ortholog in (top_orthologs[top_ortholog_seqid].keys()):
                    ortholog_gene = ".".join(ortholog.split('.')[0:-1])
                    if(plantseed_curation[ortholog_gene]['function'] == "Unannotated"):
                        continue
                    top_ortholog = ortholog
                    break
            else:
                #Ambiguously curated top orthologs, so pass
                pass

        top_ortholog_gene = ".".join(top_ortholog.split('.')[0:-1])
        if(top_ortholog_gene in plantseed_curation):
            return (top_ortholog,plantseed_curation[top_ortholog_gene]['function'],top_ortholog_seqid)
        else:
            return (None,"Unannotated",0.0)

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.token = os.environ['KB_AUTH_TOKEN']
        self.scratch = os.path.abspath(config['scratch'])
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        #END_CONSTRUCTOR
        pass

    def annotate_plant_transcripts(self, ctx, input):
        """
        :param input: instance of type "AnnotatePlantTranscriptsParams" ->
           structure: parameter "threshold" of Double, parameter "input_ws"
           of String, parameter "input_genome" of String, parameter
           "output_genome" of String
        :returns: instance of type "AnnotatePlantTranscriptsResults" ->
           structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN annotate_plant_transcripts

        output = dict()

        # Retrieve plant genome
        self.log("Fetching plant genome: "+input['input_ws']+'/'+input['input_genome'])
        plant_genome = self.dfu.get_objects({'object_refs': [input['input_ws']+'/'+input['input_genome']]})['data'][0]

        #Need to extract longest CDS, but only if CDSs available
        use_cds=1
        if(len(plant_genome['data']['cdss'])==0):
            use_cds=0
            if(len(plant_genome['data']['features'])==0):
                raise Exception("The genome does not contain any CDSs or features!")

        #Now need to be able to retrieve cds entity
        child_cds_index = dict()
        if(use_cds==1):
            child_cds_index = dict([(f['id'], i) for i, f in enumerate(plant_genome['data']['cdss'])])

        #If use_cds==1 iterate through features, iterate through CDSs, find longest sequence, use parent mRNA ID    
        #If use_cds==0 use protein_translation field if available for feature, and feature ID
        self.log("Collecting protein sequences")
        sequences_dict=dict()
        for ftr in plant_genome['data']['features']:
            if(use_cds==0 and len(ftr['protein_translation'])>0):
                sequences_dict[ftr['id']]=ftr['protein_translation']
            if(use_cds==1 and 'cdss' in ftr):
                longest_sequence=""
                longest_sequence_id=""
                for cds_id in ftr['cdss']:
                    sequence = plant_genome['data']['cdss'][child_cds_index[cds_id]]['protein_translation']
                    if(len(sequence) > len(longest_sequence)):
                        longest_sequence = sequence
                        longest_sequence_id = cds_id
                if(len(longest_sequence)>0):
                    sequences_dict[longest_sequence_id]=longest_sequence

        if(len(sequences_dict)==0):
            raise Exception("The genome does not contain any protein sequences!")

        output['ftrs'] = len(sequences_dict.keys())

        #Create directory for storing new fasta file
        uuid_string = str(uuid.uuid4())
        fasta_file_path=os.path.join(self.scratch,uuid_string)
        os.mkdir(fasta_file_path)

        #Reference data is considered immutable but each run modifies results within the directory
        #So here, we copy the reference data directory into scratch
        #The first if condition is for testing purposes
        family_file_path = ""
        testing = False
        if('families_path' in input and os.path.isdir(input['families_path'])):
            testing = True
            family_file_path = input['families_path']
            self.log("Testing Reference Families at "+family_file_path)
        else:
            uuid_string = str(uuid.uuid4())
            family_file_path=os.path.join(self.scratch,uuid_string,"Reference_Results")
            self.log("Copying Reference Families to "+family_file_path)
            shutil.copytree("/data/OrthoFinder_Phytozome_Reference",family_file_path)

        #The fasta file must have a random name to avoid _any_ clashes
        #This will need to be replaced in the newick file
        temp_genome_name = str(uuid.uuid4())
        protein_fasta_file = os.path.join(fasta_file_path,temp_genome_name+".fa")
        self.log("Printing protein sequences to file: "+protein_fasta_file)

        with open(protein_fasta_file,'w') as fasta_handle:
            #Code plagarized from https://github.com/biopython/biopython/blob/master/Bio/SeqIO/FastaIO.py
            for seq_id in sequences_dict:
                fasta_handle.write(">"+seq_id+"\n")
                for i in range(0, len(sequences_dict[seq_id]), 80):
                    fasta_handle.write(sequences_dict[seq_id][i:i+80]+"\n")

        #Building command
        command = "/kb/deployment/bin/orthofinder/orthofinder.py "
        #Software
        command +="-S diamond -M msa -A mafft -T fasttree "
        # No. of Threads
        command +="-t 8 -a 8 "
        # Avoid adding Species ID to Sequence IDs
        command += "-X "
        #For halting after alignments
        command +="-oa "
        #Input genome
        command +="-f "+fasta_file_path+" "
        #Reference families
        command +="-b "+family_file_path

        #####################################################
        output_files=list()
        if(testing is False):
            self.log("Running OrthoFinder command: "+command)

            pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

            OrthoFinder_output_file='OrthoFinder_Output.txt'
            of_fh=open(os.path.join(family_file_path,OrthoFinder_output_file),'wb')

            while pipe.poll() is None:
                stdout_line = pipe.stdout.readline()
                print(stdout_line.rstrip(), flush=True)
                of_fh.write(stdout_line)
            # Capture last piece of text if any
            stdout_line=pipe.stdout.read()
            print(stdout_line.rstrip(),flush=True)
            of_fh.write(stdout_line)
            of_fh.close()

            output_files.append({'path' : os.path.join(family_file_path,OrthoFinder_output_file),
                                 'name' : OrthoFinder_output_file,
                                 'label' : "OrthoFinder Output",
                                 'description' : 'Output text generated by OrthoFinder'})

        Ignored_Curation=dict()
        Ignored_File = "Ignored_Roles.txt"
        with open(os.path.join("/kb/module/data",Ignored_File)) as ignored_handle:
            for line in ignored_handle:
                line=line.strip()
                Ignored_Curation[line]=1

        #Parse PlantSEED families and annotation
        plantseed = FetchPlantSEEDImpl()
        plantseed_curation = plantseed.fetch_features()

        self.log("Collecting PlantSEED Curated Functions")
        PlantSEED_Roles=dict()
        output['fns']=list()
        for feature in plantseed_curation:
            function = plantseed_curation[feature]['function']
            if( function not in output['fns'] ):
                output['fns'].append(function)

            for role in plantseed_curation[feature]['roles']:
                if(role not in Ignored_Curation):
                    PlantSEED_Roles[role]=1

        output['fns']=len(output['fns'])
        output['transcripts']=len(list(plantseed_curation.keys()))
        output['alignments']=list()

        #Find, read alignments, collect families
        families_dict=dict()
        search_path = os.path.join(family_file_path,"OrthoFinder","Results_*","MultipleSequenceAlignments","OG*.fa")
        self.log("Searching for MSAs in "+search_path)
        for file in glob.glob(search_path):
            array=file.split('/')
            family=array[-1].replace(".fa","")
            alignment_handle = open(file,'rU')

            family_sequences_dict=dict()
            curated_family=list()
            # alternate header and sequence
            faiter = (x[1] for x in itertools.groupby(alignment_handle, lambda line: line[0] == ">"))
            for header in faiter:
                # drop the ">"
                header = header.__next__()[1:].strip()
                # join all sequence lines to one.
                seq = "".join(s.strip() for s in faiter.__next__())

                transcript_id=""
                try:
                    transcript_id, transcript_description = header.split(' ', 1)
                except:
                    transcript_id = header
                    transcript_description = None

                #skip un-necessary proteins to reduce computation time
                if("Athaliana" not in transcript_id and transcript_id not in sequences_dict):
                    continue

                seq = seq.upper()
                family_sequences_dict[transcript_id]=seq
                gene_id = ".".join(transcript_id.split('.')[0:-1])
                if(gene_id in plantseed_curation):
                    function = plantseed_curation[gene_id]['function']
                    curated_family.append(function+"|||"+transcript_id)
                    output['alignments'].append(transcript_id)
            
            if(len(curated_family)>0):
                if(family not in families_dict):
                    families_dict[family]=dict()
                families_dict[family]['sequences']=family_sequences_dict
                families_dict[family]['functions']=curated_family

        output['alignments']=len(output['alignments'])

        #Iterate through collected families and compute
        #Pairwise sequence identity and propagate annotation
        functions_dict=dict()
        found_annotations=list()
        clustered_features_dict=dict()
        self.log("Computing Sequence Identity on "+str(len(families_dict.keys()))+" Curated Alignments")
        for family in families_dict.keys():

            pw_seq_id_list = self.compute_clusters(families_dict[family]['sequences'])

            for function_ortholog in families_dict[family]['functions']:
                (function,ortholog) = function_ortholog.split("|||")
                if(function not in functions_dict):
                    functions_dict[function]=dict()
                if(family not in functions_dict[function]):
                    functions_dict[function][family]={'orthologs':[],
                                                      'hits':[],
                                                      'cluster':pw_seq_id_list}

                if(ortholog not in functions_dict[function][family]['orthologs']):
                    functions_dict[function][family]['orthologs'].append(ortholog)

            for spp_ftr in sorted(pw_seq_id_list.keys()):
                (ortholog,function,seqid) = self.propagate_annotation(pw_seq_id_list[spp_ftr],
                                                                      input['threshold'],
                                                                      plantseed_curation)
                ftr = spp_ftr.replace(temp_genome_name+"_","")
                clustered_features_dict[ftr]=function
                if(function not in found_annotations):
                    found_annotations.append(function)

                if(function in functions_dict and \
                       family in functions_dict[function] and \
                       ortholog in functions_dict[function][family]['orthologs']):
                    functions_dict[function][family]['hits'].append({'seqid':seqid,
                                                                     'feature':spp_ftr})
        output['hit_fns']=len(found_annotations)
        output['hit_ftrs']=len(clustered_features_dict.keys())

        #Now, re-populate feature functions, and save genome object
        #But, if annotating CDS, need to be able to retrieve parent feature
        parent_feature_index = dict()
        if(use_cds==1):
            parent_feature_index = dict([(f['id'], i) for i, f in enumerate(plant_genome['data']['features'])])

        self.log("Populating plant genome with newly clustered functions")
        #Add annotation to protein-coding genes
        for ftr in plant_genome['data']['features']:
            ftr['functions']=["Unannotated"]
            if(ftr['id'] in clustered_features_dict):
                ftr['functions']=[clustered_features_dict[ftr['id']]]

        #Add annotation to transcripts
        for mrna in plant_genome['data']['mrnas']:
            mrna['functions']=["Unannotated"]
            if(mrna['id'] in clustered_features_dict):
                mrna['functions']=[clustered_features_dict[mrna['id']]]
                plant_genome['data']['features'][parent_feature_index[mrna['parent_gene']]]['functions']=[clustered_features_dict[mrna['id']]]
            elif(mrna['parent_gene'] in clustered_features_dict):
                mrna['functions']=[clustered_features_dict[mrna['parent_gene']]]
                plant_genome['data']['features'][parent_feature_index[mrna['parent_gene']]]['functions']=[clustered_features_dict[mrna['parent_gene']]]

        #Add annotation to proteins
        for cds in plant_genome['data']['cdss']:
            cds['functions']=["Unannotated"]
            if(cds['id'] in clustered_features_dict):
                cds['functions']=[clustered_features_dict[cds['id']]]
                plant_genome['data']['features'][parent_feature_index[cds['parent_gene']]]['functions']=[clustered_features_dict[cds['id']]]
            elif(cds['parent_gene'] in clustered_features_dict):
                cds['functions']=[clustered_features_dict[cds['parent_gene']]]
                plant_genome['data']['features'][parent_feature_index[cds['parent_gene']]]['functions']=[clustered_features_dict[cds['parent_gene']]]
                
        #Save genome
        with open("/kb/module/work/tmp/shit.json","w") as fh:
            import json
            json.dump(plant_genome, fh)
            fh.close()

        if('output_genome' not in input):
            input['output_genome']=input['input_genome']

        #wsid = self.dfu.ws_name_to_id(input['input_ws'])
        #save_result = self.dfu.save_objects({'id':wsid,'objects':[{'name':input['output_genome'],
        #                                                           'data':plant_genome['data'],
        #                                                           'type':'KBaseGenomes.Genome'}]})[0]
        
        save_result = self.gfu.save_one_genome({'workspace' : input['input_ws'],
                                                'name' : input['output_genome'],
                                                'data' : plant_genome['data']})['info'];

        #reference of saved genome
        saved_genome = "{}/{}/{}".format(save_result[6],save_result[0],save_result[4])

        Annotated_Roles=dict()
        for curation in found_annotations:
            Function_Comments = curation.split("#")
            for i in range(len(Function_Comments)):
                Function_Comments[i]=Function_Comments[i].strip()

            Function = Function_Comments.pop(0)
            Roles = re.split("\s*;\s+|\s+[\@\/]\s+", Function)
            for role in Roles:
                if(role in PlantSEED_Roles):
                    Annotated_Roles[role]=1

        output['hit_fns']=len(found_annotations)
        output['cur_roles']=len(PlantSEED_Roles.keys())
        output['hit_roles']=len(Annotated_Roles.keys())

        # Calculate fraction of PlantSEED functional roles
        fraction_plantseed = float( (float(len(Annotated_Roles.keys())) / float(len(PlantSEED_Roles.keys()))) )

        # HTML Folder Path
        uuid_string = str(uuid.uuid4())
        html_file_path=os.path.join(self.scratch,uuid_string)
        os.mkdir(html_file_path)

        # Generate figure: 
        #     the path parameter is for the reference data that is integrated into the figure.
        figure_generator = GenerateFigureImpl("/kb/module/data/")
        bokeh_figure = figure_generator.generate_figure(input['threshold'],
                                                        fraction_plantseed)

        # Save figure
        figure_html_file="figure.html"
        output_file(os.path.join(html_file_path,figure_html_file))
        save(bokeh_figure)

        # Generate table
        table_generator = GenerateTableImpl()
        annotation_table_string = table_generator.generate_table(functions_dict)

        # Save table
        # Read in template html
        with open(os.path.join('/kb/module/data',
                               'app_report_templates',
                               'annotation_report_tables_template.html')) as report_template_file:
            report_template_string = report_template_file.read()

        # Generate and Insert html title
        #     This needs to be done because it affects the name of the CSV download
        title_string = "-".join([input['input_genome'],str(input['threshold'])])
        report_template_string = report_template_string.replace('*TITLE*', title_string)

        # Insert table into template
        table_report_string = report_template_string.replace('*TABLES*', annotation_table_string)

        # Save table
        table_html_file = "table.html"
        with open(os.path.join(html_file_path,table_html_file),'w') as table_file:
            table_file.write(table_report_string)

        # Generate main index.html content
        html_string="<html><head><title>KBase Plant OrthoFinder Report</title></head><body>"
        html_string+="<div style=\"text-align: center; max-width: 800px\">"
        html_string+="<p>The Plant OrthoFinder app has finished running: "
        html_string+=str(output['ftrs'])+" protein sequences were clustered "
        html_string+="with "+str(output['transcripts'])+ " PlantSEED-curated enzymes. "
        html_string+="The app was able to predict "+str(output['hit_fns'])+" enzymatic functions "
        html_string+="for "+str(output['hit_ftrs'])+" protein sequences and "
        html_string+="this result indicates that, for this set of protein sequences, "
        html_string+="the app detected {0:.0f}%".format(float(fraction_plantseed*100.0))
        html_string+=" of the enzymatic functions of plant primary metabolism that were "
        html_string+="curated as part of the PlantSEED project.</p></br>"
        html_string+="<p>The results of the annotation are tabulated in this "
        html_string+="<a href=\""+table_html_file+"\" target=\"_blank\">Table</a></p></div>"

        caption="<figcaption><b>Figure 1: Propagation of metabolic roles for "
        caption+=str(input['input_genome'])+". </b>"
        caption+="The PlantSEED project curated "+str(output['cur_roles'])
        caption+=" distinct primary metabolic roles for Arabidopsis thaliana. "
        caption+="Here we show the impact of propagating these roles to other "
        caption+="species using sequence identity. "
        caption+="For each group of species, and for a different threshold of "
        caption+="sequence identity, we show the fraction of curated roles that were "
        caption+="propagated. The fraction of propagated roles decreases as a function "
        caption+="of similarity and phylogenetic distance. "
        caption+="The fraction of roles that were propagated for "+input['input_genome']
        caption+=" at the chosen threshold of "
        caption+=str(input['threshold'])+" for sequence identity is {0:.2f}".format(float(fraction_plantseed))
        caption+=" and is marked by the bold plus. "
        caption+=" A user may re-run the app with a different threshold; a higher "
        caption+="threshold will increase the reliability of the results, but will "
        caption+=" reduce the number of propagated metabolic roles; a lower threshold "
        caption+="will increase the number of propagated metabolic roles but "
        caption+=" also increase the number of propagations, increasing the likelihood "
        caption+="of a false positive. This figure can be viewed in a separate window "
        caption+="<a href=\""+figure_html_file+"\" target=\"_blank\">here</a></figcaption>"

        html_string+="<div style=\"text-align: center; max-width: 620px\">"
        html_string+="<embed type=\"text/html\" src=\""+figure_html_file+"\" width=\"620\" height=\"620\"></embed>"
        html_string+=caption
        html_string+="</div></body></html>"

        #Save index file
        with open(os.path.join(html_file_path,"index.html"),'w') as index_file:
            index_file.write(html_string)

        upload_info = self.dfu.file_to_shock({'file_path': html_file_path,
                                              'pack': 'zip'})

        html_report_list=list()
        html_link = {'shock_id' : upload_info['shock_id'],
                     'name' : figure_html_file,
                     'label' : 'Scatterplot figures generated by Integrate Abundances with Metabolism app',
                     'description' : 'Scatterplot figures generated by Integrate Abundances with Metabolism app'}
        html_report_list.append(html_link)

        html_link = {'shock_id' : upload_info['shock_id'],
                     'name' : table_html_file,
                     'label' : 'Scatterplot figures generated by Integrate Abundances with Metabolism app',
                     'description' : 'Scatterplot figures generated by Integrate Abundances with Metabolism app'}
        html_report_list.append(html_link)

        description = "Plant genome "+plant_genome['data']['id']+" annotated with metabolic functions"

        uuid_string = str(uuid.uuid4())
        report_params = { 'direct_html_link_index' : 0, #Use to refer to index of 'html_links'
                          'workspace_name' : input['input_ws'],
                          'report_object_name' : 'kb_orthofinder_' + uuid_string,
                          'objects_created' : [{"ref":saved_genome,"description":description}],
                          'file_links' : output_files,
                          'html_links' : html_report_list}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        report_client_output = kbase_report_client.create_extended_report(report_params)
        output['report_name']=report_client_output['name']
        output['report_ref']=report_client_output['ref']

        #END annotate_plant_transcripts

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method annotate_plant_transcripts return value ' +
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
