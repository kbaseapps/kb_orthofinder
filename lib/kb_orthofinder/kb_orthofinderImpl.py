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

from KBaseReport.KBaseReportClient import KBaseReport
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as workspaceService
from kb_orthofinder.GenerateFigure import GenerateFigure

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
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:kbaseapps/kb_orthofinder.git"
    GIT_COMMIT_HASH = "3ca54e22b7804ac0582589aa185e54fcdf0b48b1"

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

                Spp1=features[i].split("|_|")[0]
                Spp2=features[j].split("|_|")[0]
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

    def propagate_annotation(self, family, cluster, threshold, arabidopsis_functions):
        top_orthologs = dict()
        for ortholog in (cluster.keys()):
            
            if(cluster[ortholog] not in top_orthologs):
                top_orthologs[cluster[ortholog]]=dict()
            top_orthologs[cluster[ortholog]][ortholog]=1

        top_ortholog_seqid=0.0
        for seq_id in (top_orthologs.keys()):
            if(seq_id>top_ortholog_seqid):
                top_ortholog_seqid=seq_id

        #Key to controlling the propagation of annotation
        #Is the sequence identity of the highest-scoring ortholog good enough?
        if(top_ortholog_seqid < threshold):
            return "Uncurated"

        top_ortholog=""
        if(len(top_orthologs[top_ortholog_seqid])==1):
            top_ortholog = top_orthologs[top_ortholog_seqid].keys()[0]
        else:
            #There are multiple Arabidopsis orthologs that have the same level of seq. id
            #Collect the actual functions and see if there's multiple functions

            Multi_Functions = dict()
            for ortholog in (top_orthologs[top_ortholog_seqid].keys()):
                if(ortholog not in arabidopsis_functions):
                    print family,ortholog
                    function = "Uncurated"
                    arabidopsis_functions[ortholog]=function
                else:
                    function = arabidopsis_functions[ortholog]
                Multi_Functions[function]=1

            #Rules are:
            # 0) if 1 function, arbitrarily pick ath ortholog
            # 1) if 2 functions, and one is "Uncurated", prioritize curation (arbitrarily pick ath ortholog that is curated)
            # 2) if 2 functions, and none is "Uncurated" or > 2 functions, make it ambiguous

            if(len(Multi_Functions.keys())==1):
                top_ortholog = top_orthologs[top_ortholog_seqid].keys()[0]
            elif(len(Multi_Functions.keys())==2 and "Uncurated" in Multi_Functions.keys()):
                for ortholog in (top_orthologs[top_ortholog_seqid].keys()):
                    if(arabidopsis_functions[ortholog] == "Uncurated"):
                        continue
                    top_ortholog = ortholog
                    break
            else:
                #Ambiguously curated top orthologs, so pass
                pass

        if(top_ortholog in arabidopsis_functions):
            return arabidopsis_functions[top_ortholog]
        else:
            return "Uncurated"

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

        # Force upgrade
        if("feature_counts" in plant_genome['data']):
            del(plant_genome['data']['feature_counts'])

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
            if(use_cds==1):
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
        self.log("Copying Reference Families")
        family_file_path = ""
        if('families_path' in input and os.path.isdir(input['families_path'])):
            family_file_path = input['families_path']
        else:
            uuid_string = str(uuid.uuid4())
            family_file_path=os.path.join(self.scratch,uuid_string,"Reference_Results")
            shutil.copytree("/data/Reference_Results",family_file_path)

        #File handle
        self.log("Printing protein sequences to file")
        #The fasta file must have a random name to avoid _any_ clashes
        #This will need to be replaced in the newick file
        temp_genome_name = str(uuid.uuid4())
        with open(os.path.join(fasta_file_path,temp_genome_name+".fa"),'w') as fasta_handle:
            #Code plagarized from https://github.com/biopython/biopython/blob/master/Bio/SeqIO/FastaIO.py
            for seq_id in sequences_dict:
                fasta_handle.write(">"+seq_id+"\n")
                for i in range(0, len(sequences_dict[seq_id]), 80):
                    fasta_handle.write(sequences_dict[seq_id][i:i+80]+"\n")

        #Building command
        command = "/kb/deployment/bin/orthofinder/orthofinder.py "
        #Software
        command +="-S diamond -M msa -A mafft -T fasttree "
        #Threads
        command +="-t 8 -a 8 "
        #For halting after alignments
        command +="-oa "
        #Input genome
        command +="-f "+fasta_file_path+" "
        #Reference families
        command +="-b "+family_file_path

        self.log("Running OrthoFinder command: "+command)
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

        OrthoFinder_output_file='OrthoFinder_Output.txt'
        of_fh=open(os.path.join(family_file_path,OrthoFinder_output_file),'w')

        while pipe.poll() is None:
            stdout_line = pipe.stdout.readline()
            print stdout_line.rstrip()
            of_fh.write(stdout_line)
        #Capture last piece of text if any
        stdout_line=pipe.stdout.read()
        print stdout_line.rstrip()
        of_fh.write(stdout_line)
        of_fh.close()

#        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#        output_tuple = pipe.communicate()
#        exitCode = pipe.returncode

        output_files=list()
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
        PlantSEED_Curation=dict()
        Curation_File = "Arabidopsis_Family_Curation.txt"

        output['fns']=dict()
        self.log("Collecting PlantSEED Curation")
        PlantSEED_Roles=dict()
        with open(os.path.join("/kb/module/data",Curation_File)) as plantseed_families_handle:
            for line in plantseed_families_handle.readlines():
                line=line.strip()
                (transcript,curation)=line.split('\t')
                PlantSEED_Curation[transcript]=curation
                output['fns'][curation]=1
         
                #Split out comments
                Function_Comments = curation.split("#")
                for i in range(len(Function_Comments)):
                    Function_Comments[i]=Function_Comments[i].strip()

                Function = Function_Comments.pop(0)
                Roles = re.split("\s*;\s+|\s+[\@\/]\s+", Function)
                for role in Roles:
                    if(role not in Ignored_Curation):
                        PlantSEED_Roles[role]=1

        output['fns']=len(output['fns'].keys())
        output['transcripts']=len(PlantSEED_Curation.keys())
        output['alignments']=list()

        #Find, read alignments, collect families
        families_dict=dict()
        self.log("Searching for MSAs")
        for file in glob.glob(os.path.join(family_file_path,"Orthologues_*","Alignments","OG*.fa")):
            array=file.split('/')
            family=array[-1].replace(".fa","")
            alignment_handle = open(file,'rU')

            family_sequences_dict=dict()
            curated_family=0
            # alternate header and sequence
            faiter = (x[1] for x in itertools.groupby(alignment_handle, lambda line: line[0] == ">"))
            for header in faiter:
                # drop the ">"
                header = header.next()[1:].strip()
                # join all sequence lines to one.
                seq = "".join(s.strip() for s in faiter.next())

                fasta_header=""
                try:
                    fasta_header, fasta_description = header.split(' ', 1)
                except:
                    fasta_header = header
                    fasta_description = None

                #skip un-necessary proteins to reduce computation time
                spp=fasta_header.split("|_|")[0]
                if("Athaliana" not in spp and temp_genome_name not in spp):
                    continue

                #This should've been un-necessary but here I need to remove a redundancy
                fasta_header = fasta_header.replace("Athaliana_TAIR10_Athaliana","Athaliana")
                
                seq = seq.upper()
                family_sequences_dict[fasta_header]=seq
                if(fasta_header in PlantSEED_Curation):
                    curated_family=1
                    output['alignments'].append(fasta_header)
            
            if(curated_family==1):
                families_dict[family]=family_sequences_dict

        output['alignments']=len(output['alignments'])

        #Iterate through collected families and compute
        #Pairwise sequence identity and propagate annotation
        found_annotations_dict=dict()
        clustered_features_dict=dict()
        self.log("Computing Sequence Identity on "+str(len(families_dict.keys()))+" Curated Alignments")
        for family in families_dict.keys():
            pw_seq_id_list = self.compute_clusters(families_dict[family])
            for spp_ftr in sorted(pw_seq_id_list.keys()):
                annotation = self.propagate_annotation(family,
                                                       pw_seq_id_list[spp_ftr],
                                                       input['threshold'],
                                                       PlantSEED_Curation)
                ftr = spp_ftr.replace(temp_genome_name+"_","")
                clustered_features_dict[ftr]=annotation
                found_annotations_dict[annotation]=1

        output['hit_fns']=len(found_annotations_dict.keys())
        output['hit_ftrs']=len(clustered_features_dict.keys())

        #Now, re-populate feature functions, and save genome object
        #But, if annotating CDS, need to be able to retrieve parent feature
        parent_feature_index = dict()
        if(use_cds==1):
            parent_feature_index = dict([(f['id'], i) for i, f in enumerate(plant_genome['data']['features'])])

        self.log("Populating plant genome with newly clustered functions")
        #Add annotation to protein-coding genes
        for ftr in plant_genome['data']['features']:
            ftr['functions']=["Uncurated"]
            if(ftr['id'] in clustered_features_dict):
                ftr['functions']=[clustered_features_dict[ftr['id']]]

        #Add annotation to proteins
        for cds in plant_genome['data']['cdss']:
            cds['functions']=["Uncurated"]
            if(cds['id'] in clustered_features_dict):
                cds['functions']=[clustered_features_dict[cds['id']]]
                plant_genome['data']['features'][parent_feature_index[cds['parent_gene']]]['functions']=[clustered_features_dict[cds['id']]]
                
        #Save genome
        if('output_genome' not in input):
            input['output_genome']=input['input_genome']

        save_result = self.gfu.save_one_genome({'workspace' : input['input_ws'],
                                                'name' : input['output_genome'],
                                                'data' : plant_genome['data']});

        #Calculate fraction of PlantSEED functional roles
        Annotated_Roles=dict()
        for curation in found_annotations_dict.keys():
            #Split out comments
            Function_Comments = curation.split("#")
            for i in range(len(Function_Comments)):
                Function_Comments[i]=Function_Comments[i].strip()

            Function = Function_Comments.pop(0)
            Roles = re.split("\s*;\s+|\s+[\@\/]\s+", Function)
            for role in Roles:
                if(role in PlantSEED_Roles):
                    Annotated_Roles[role]=1

        output['cur_roles']=len(PlantSEED_Roles.keys())
        output['hit_roles']=len(Annotated_Roles.keys())
        fraction_plantseed = float( (float(len(Annotated_Roles.keys())) / float(len(PlantSEED_Roles.keys()))) )
        figure_params = {'threshold':input['threshold'],'fraction':fraction_plantseed,'id':input['input_genome']}
        figure_path = self.generate_figure(figure_params)

        html_string="<html><head><title>KBase Plant OrthoFinder Report</title></head><body>"
        html_string+="<p>The Plant OrthoFinder app has finished running. "
        html_string+=str(output['ftrs'])+" protein sequences were clustered with "+str(output['transcripts'])+ " PlantSEED-curated enzymes.</p>"
        html_string+="<p>The app was able to predict "+str(output['hit_fns'])+" enzymatic functions for "+str(output['hit_ftrs'])+" protein sequences.</p>"
        html_string+="<p>This result indicates that, for this set of protein sequences, the app detected {0:.0f}%".format(float(fraction_plantseed*100.0))
        html_string+=" of the enzymatic functions of plant primary metabolism that were curated as part of the PlantSEED project.</p>"

        caption="<figcaption><b>Figure 1: </b>The PlantSEED project curated str(output['roles']) distinct enzymatic roles for Arabidopsis thaliana. "
        caption+="Here we show the impact of propagating these enzymatic roles to other species using sequence identity. "
        caption+="For each group of species, and for a different threshold of sequence identity, we show the fraction of curated roles that were "
        caption+="propagated. The fraction of roles that were propagated for "+input['input_genome']+" at the chosen threshold of "
        caption+=str(input['threshold'])+" for sequence identity is {0:.0f}%".format(float(fraction_plantseed))
        caption+=" and is marked by the bold plus.</figcaption>"

        for file in os.listdir(figure_path):
            format = file.split('.')[-1].upper()

            if(format == "PNG"):
#                html_string+="<p><img src=\""+file+"\"/></p></body></html>"
                html_string+="<center><figure><img src=\""+file+"\"/>"+caption+"</figure></center>"

            output_files.append({'path' : os.path.join(figure_path,file),
                                 'name' : file,
                                 'label' : "PyGrace Figure",
                                 'description' : 'PyGrace Figure in '+format+' format'})

        #Save index file
        with open(os.path.join(figure_path,"index.html"),'w') as index_file:
            index_file.write(html_string)

        upload_info = self.dfu.file_to_shock({'file_path': figure_path,
                                              'pack': 'zip'})
            
        html_folder = {'shock_id' : upload_info['shock_id'], #Used instead of 'path'
                       #'path' : figure_path, #KBaseReport zip_archive() is broken
                       'name' : 'index.html', #Used for URL path
                       'label' : 'html files',
                       'description' : 'HTML files'}

        saved_genome = "{}/{}/{}".format(save_result['info'][6],save_result['info'][0],save_result['info'][4])
        description = "Plant genome "+plant_genome['data']['id']+" annotated with metabolic functions"

        uuid_string = str(uuid.uuid4())
        report_params = { 'objects_created' : [{"ref":saved_genome,"description":description}],
                          'file_links' : output_files,
                          'html_links' : [html_folder],
                          'direct_html_link_index' : 0, #Use to refer to index of 'html_links'
#                          'direct_html' : html_string, # Can't embed images
                          'workspace_name' : input['input_ws'],
                          'report_object_name' : 'kb_plant_rast_report_' + uuid_string }
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
