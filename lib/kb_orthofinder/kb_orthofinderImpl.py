# -*- coding: utf-8 -*-
#BEGIN_HEADER
import glob
import os
import uuid
import shutil
import subprocess
import itertools

from KBaseReport.KBaseReportClient import KBaseReport
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as workspaceService
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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER

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
        family_file_path = ""
        if(input['families_path']):
            family_file_path = input['families_path']
        else:
            uuid_string = str(uuid.uuid4())
            family_file_path=os.path.join(self.scratch,uuid_string,"Reference_Results")
            shutil.copytree("/data/Reference_Results",family_file_path)

        #File handle
        with open(os.path.join(fasta_file_path,input['input_genome']+'.fa'),'w') as fasta_handle:
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
        #Input genome
        command +="-f "+fasta_file_path+" "
        #Reference families
        command +="-b "+family_file_path

        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) 
        output_string = pipe.communicate()[0] 
        exitCode = pipe.returncode

        if (exitCode != 0): 
            error_msg = "Error running command: \n"+command+"\n"
            error_msg +="Exit Code: "+str(exitCode)+"\n"
            error_msg +="Output: \n"+output+"\n"
            raise ValueError(error_msg)

        report_string = "Executed command: \n"+command+"\n"
        report_string +="Exit Code: "+str(exitCode)+"\n"
        report_string +="Output: \n"+output_string+"\n"
        output['report']=report_string

        #Parse PlantSEED families and annotation
        PlantSEED_Curation=dict()
        Curation_File = "Arabidopsis_Family_Curation.txt"

        output['fns']=dict()
        with open(os.path.join("/kb/module/data",Curation_File)) as plantseed_families_handle:
            for line in plantseed_families_handle.readlines():
                line=line.strip()
                (transcript,curation)=line.split('\t')
                PlantSEED_Curation[transcript]=curation
                output['fns'][curation]=1
         
        output['transcripts']=PlantSEED_Curation.keys()
        output['alignments']=list()

        #Find, read alignments, collect families
        families_dict=dict()
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
                if("Athaliana" not in spp and input['input_genome'] not in spp):
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

        #Iterate through collected families and compute
        #Pairwise sequence identity and propagate annotation
        found_annotations_dict=dict()
        clustered_features_dict=dict()
        for family in families_dict.keys():
            pw_seq_id_list = self.compute_clusters(families_dict[family])
            for spp_ftr in sorted(pw_seq_id_list.keys()):
                annotation = self.propagate_annotation(family,
                                                       pw_seq_id_list[spp_ftr],
                                                       input['threshold'],
                                                       PlantSEED_Curation)
                ftr = spp_ftr.replace(input['input_genome']+"_","")
                clustered_features_dict[ftr]=annotation
                found_annotations_dict[annotation]=1

        output['hit_fns']=len(found_annotations_dict.keys())
        output['hit_ftrs']=len(clustered_features_dict.keys())

        #Now, re-populate feature functions, and save genome object
        #But, if annotating CDS, need to be able to retrieve parent feature
        parent_feature_index = dict()
        if(use_cds==1):
            parent_feature_index = dict([(f['id'], i) for i, f in enumerate(plant_genome['data']['features'])])

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

        html_string="<html><head><title>KBase Plant OrthoFinder Report</title></head><body>"
        html_string+="<p>The Plant OrthoFinder app has finished running. "
        html_string+=str(output['ftrs'])+" protein sequences were clustered with "+str(len(output['transcripts']))+ " PlantSEED-curated enzymes.</p>"
        html_string+="<p>The app was able to predict "+str(output['hit_fns'])+" enzymatic functions for "+str(output['hit_ftrs'])+" protein sequences.</p>"
        fraction_plantseed = float( (float(output['hit_fns']) / float(len(output['fns'].keys()))) * 100.0 )
        html_string+="<p>This result indicates that, for this set of protein sequences, the app detected {0:.0f}%".format(fraction_plantseed)
        html_string+=" of the enzymatic functions of plant primary metabolism that were curated as part of the PlantSEED project.</p></body></html>"

        print html_string

        saved_genome = "{}/{}/{}".format(save_result['info'][6],save_result['info'][0],save_result['info'][4])
        description = "Plant genome "+plant_genome['data']['id']+" annotated with metabolic functions"
        uuid_string = str(uuid.uuid4())
        report_params = { 'objects_created' : \
                          [{"ref":saved_genome,"description":description}],
                          'direct_html' : html_string,
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
