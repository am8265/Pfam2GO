#!/usr/bin/env python
import re
import sys
import math
def gid_pfam(fh_1,eval=math.exp(-10)):
    '''gid_pfam(hmmscan output in --domtblout format,default evalue is 1e-10)\
    the function returns a dictionary of dictinaries: gid_pfam_evalues={gid:{(pfam1,pfam1_des):[evals],(pfam2,pfam2_des):[evals]...(pfamn,pfamn_des):[evals]}'''
    fh=open(fh_1)
    gid_pfam_eval=dict()
    for line in fh:
        line=line.rstrip('\n').strip()
        if re.search('(PF\d+)\.\d+',line)!=None:
            pfam_gid=re.search('(PF\d+)\.\d+\s+(\d+)\s*(\S+)\s+\S+\s+\S+\s+(\S+).+\d{0,1}\.\d{0,3}\S(.+)',line).groups()#This is where we pick up 


            pfam_id=pfam_gid[0]
        #print(pfam_id)
            gid=pfam_gid[2]
        #print(gid)
            def_pfam=pfam_gid[4]#pfam_definition
            evalue=float(pfam_gid[3])
######a single UNIQUE gid or peptide_ID will have (pfam_id,pfam_def, evalue)....It could potentially have SAME (pfam_id,pfam_def withe different evalue)
            if evalue<=eval:#if there is a cutoff match....if the evalue cut off is below the listed threshold
                if gid_pfam_eval.get(gid,'0')=='0':#i.e we are encountering the PEPTIDE_id 1st time!
                    gid_pfam_eval[gid]={(pfam_id,def_pfam):[evalue]}
                elif gid_pfam_eval[gid].get((pfam_id,def_pfam),'NA')=='NA':
                    gid_pfam_eval[gid][(pfam_id,def_pfam)]=[evalue]
                else:
                    gid_pfam_eval[gid][(pfam_id,def_pfam)].append(evalue)#we append if the gid already exist!---Remember a given gid could have several same or different  PFAM domains distributed
                    gid_pfam_eval[gid][(pfam_id,def_pfam)]=list(set(gid_pfam_eval[gid][(pfam_id,def_pfam)]))
            else:
                pass#we we just ignore the line if its above the particular evalue threshold 
    fh.close() 
    return gid_pfam_eval       

###########
out_fh=open('gid_pfam_def.txt','w')
out_fh.write('peptide_id'+'\t'+'pfam_id'+'\t'+'definition'+'\t'+'Evalues'+'\n')
gid_pfam_eval=gid_pfam(sys.argv[1])
######running diagnistics!!!!!!
print(gid_pfam_eval)
########writing to the gid_pfam_def.txt file####################################
for ids in sorted(gid_pfam_eval.keys()):
    pfam_id_cat=''
    pfam_def_cat=''
    eval_min_cat=''
    for pfam,evals in gid_pfam_eval[ids].items():
        pfam_id=pfam[0].strip().strip("'")
        pfam_def=pfam[1].strip().strip("'")
        eval_min=str(min(evals))
        pfam_id_cat+=pfam_id+','
        pfam_def_cat+=pfam_id+'['+pfam_def+']; '
        eval_min_cat+=pfam_id+'['+eval_min+']; '
    pfam_id_cat=pfam_id_cat.strip().rstrip(',')
    pfam_def_cat=pfam_def_cat.strip().rstrip(';')
    eval_min_cat=eval_min_cat.strip().rstrip(';')
    out_fh.write(ids+'\t'+pfam_id_cat+'\t'+pfam_def_cat+'\t'+eval_min_cat+'\n')
gid_pfam_eval=gid_pfam(sys.argv[1])
print(gid_pfam_eval)

