#!/usr/bin/env python 

import re
import os
import math
import subprocess


class Pfam2GO:

    def run_pfam (self, pfam_dir, protein_fasta_file):
    	
        """
        _run_pfam: scan the protein_fasta_file for PFAM domain using pfam-A.hmm INDEXED files (.hmm files pressed using)present in the
        pfam_dir(pfam_dir = /<path_to>/Pfam-A.hmm) and returns a file name hmmscan_file in the working directory
        """
        #define path to hmmscan binaries
        hmmscan_path = '/Users/ayanmalakar/BioinformaticsCookBook/Pfam_GO/gid_go_def/hmmer-3.1b2-macosx-intel/binaries/hmmscan'
        completed = subprocess.run([hmmscan_path,'--domtblout','hmmscan_file',pfam_dir,protein_fasta_file], check="TRUE", stdout=subprocess.PIPE)

        return './hmmscan_file'#defined output file path



    	
    def map_geneid_pfam(self,hmmscan_file,e=math.exp(-10)):

        gid_pfam_eval=dict()
        gid_pfam_def=dict()
        hmmscan_file=open(hmmscan_file)
        for line in hmmscan_file:

            line=line.rstrip('\n').strip()
            if re.search('(PF\d+)\.\d+',line)!=None:

                pfam_gid=re.search('(PF\d+)\.\d+\s+(\d+)\s*(\S+)\s+\S+\s+\S+\s+(\S+).+\d{0,1}\.\d{0,3}\S(.+)',line).groups()

        #PFAM ID

                pfam_id=pfam_gid[0]
                #Gene Identifier

                geneid=pfam_gid[2]
                #PFAM DEFINITION

                pfam_def=pfam_gid[4]

                #E-VALUE for the PFAM DOMAIN IN "HMMSCAN FILE"
                evalue=float(pfam_gid[3])

                if evalue <= e: #Accept maximum E-values for PFAM domain hits  to that of user provided 'e' or E-value

                    if gid_pfam_eval.get(geneid,'0')=='0': # i.e we are encountering the Peptide_id/GeneID 1st time!

                        gid_pfam_eval[geneid]=[(pfam_id, evalue)]

                        gid_pfam_def[geneid]=[(pfam_id, pfam_def)]

                    else :

                        gid_pfam_eval[geneid].append((pfam_id, evalue))
    #remove any duplicate (PFAM,eval) pair
                        gid_pfam_eval[geneid]=list(set(gid_pfam_eval[geneid]))


                        gid_pfam_def[geneid].append((pfam_id, pfam_def))
    #remove any duplicate (PFAM,PFAM_definition) pair
                        gid_pfam_def[geneid]=list(set(gid_pfam_def[geneid]))

                else :
        #we we just ignore the entries if its above the particular e-value threshold
                    pass

        hmmscan_file.close()


        return (gid_pfam_eval,gid_pfam_def)

    


    def pfam_2_go(self,pfamtoGO_mapping_file):

        pfam_2_GO=dict()
        pfamtoGO_mapping_file=open(pfamtoGO_mapping_file)
        for line in pfamtoGO_mapping_file:
            #print(line)
            if re.search('^Pfam',line)!=None:

                line=line.rstrip('\n').strip()
                parse_line=line.split(':')[1].split(' ')
                pfam=parse_line[0]
            #print(pfam)
                parse_line=':'.join(line.split(':')[2:])
                GO=parse_line.split(';')[1].strip().strip(' ')
                #print(pfam)
                if pfam_2_GO.get(pfam,0)==0:#if PFAM_id doesn't exist then create a new one!
                #pfam2go[pfam]=go
                    pfam_2_GO[pfam]=[GO]
            #print(pfam2go)

                else:#if the PFAM id already exist as a key then append GO Terms to previous list
    #pfam2go[pfam]=pfam2go[pfam]+','+go
                    pfam_2_GO[pfam].append(GO)
                #Removing the possibility of a DUPLICATE GO
                    pfam_2_GO[pfam]=list(set(pfam_2_GO[pfam]))

        pfamtoGO_mapping_file.close()

        return pfam_2_GO



    def gene_2_go (self,gid_pfam_eval,pfam_2_GO):
        geneid_2_GO=dict()
        for gid,pfam_eval in gid_pfam_eval.items():
            pfam_list = list(set([pfam[0] for pfam in pfam_eval]))
            geneid_2_GO[gid]=[]
            for pfam in pfam_list:
                if pfam_2_GO.get(pfam, "NA")!='NA': #i.e we have a pfam with GOs info
                    geneid_2_GO[gid].extend(pfam_2_GO[pfam])
                else:#We have pfams without corresponding GOs then Don't do anything
                    pass
            geneid_2_GO[gid]=list(set(geneid_2_GO[gid])) # removed OVERLAPPING GO terms for the geneID(gid)
        
        #Its possible geneid_2_GO[gid]=[] for cases where None of the corresponding PFAM had ANY defined GO from pfam_2_GO
        return geneid_2_GO

    def write_2_files(self,gid_pfam_def,geneid_2_GO):
        fh1 = open("gid2PFAM.tsv", 'w')
        fh1.write("feature_ID"+'\t'+'PFAM_ID'+'\t'+'PFAM_Def'+'\n')
        for gid,pfam_def in gid_pfam_def.items():
            for pfam,define in pfam_def:
                fh1.write(gid+'\t'+pfam+'\t'+define+'\n')
        fh1.close()

        fh2=open("gid2GO.tsv", 'w')
        fh2.write("feature_ID"+'\t'+'GO_ID'+'\n')
        for gid,GOs in geneid_2_GO.items():
            for GO in GOs:
                print(GO)
                fh2.write(gid+'\t'+GO+'\n')
        fh2.close()


    def __init__(self,config):

        pass


    def run_pfam_2_go(self,pfam_dir,protein_fasta_file,pfamtoGO_mapping_file):
        hmmscan_file=self.run_pfam(pfam_dir,protein_fasta_file)
        gid_pfam_eval,gid_pfam_def=self.map_geneid_pfam(hmmscan_file)#evalue_cutoff
        pfam_2_GO=self.pfam_2_go(pfamtoGO_mapping_file)
        print(pfam_2_GO)
        gene_2_GO=self.gene_2_go(gid_pfam_eval,pfam_2_GO)
        print(gene_2_GO)
#we write to predefined files ---could possibly be user provided path 
        self.write_2_files(gid_pfam_def,gene_2_GO)
        return gid_pfam_def,gene_2_GO








