#!/usr/bin/env python
import re
import sys
import math
def gid_pfam(fh_1,eval=math.exp(-10)):
    
    fh=open(fh_1)
    #gid_pfam=dict()
    #pfam_def=dict()
    gid_pfam_eval=dict()
    for line in fh:
        line=line.rstrip('\n').strip()
        if re.search('(PF\d+)\.\d+',line)!=None:

######## Rewrite the parser script##################
########Become knowledgable in github############
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
#     for ids in sorted(gid_pfam_eval.keys()):
         

            
#                if pfam_def.get(pfam_id,'0')=='0':
#                    pfam_def[pfam_id]=[def_pfam]#every unique pfam_id shoudl have a unique definition
#                else:
#                    pfam_def[pfam_id].append(def_pfam)
            
#                if gid_pfam.get(gid,'0')=='0':
#                    gid_pfam[gid]=[pfam_id]
#                else:
#                    gid_pfam[gid].append(pfam_id)
    
#print(gid_pfam)
#    for ids in sorted(gid_pfam.keys()):
#        gid_pfam[ids]=list(set(gid_pfam[ids]))
#    for ids in pfam_def.keys():
#        pfam_def[ids]=list(set(pfam_def[ids]))
#    return (gid_pfam,pfam_def)
#    fh.close()

gid_pfam_eval=gid_pfam(sys.argv[1])
print(gid_pfam_eval)
 




#id_pfam=gid_pfam[0]
#pfam_def=gid_pfam[1]
######testing here
#print(id_pfam)
##########
out_fh=open('gid_pfam_def.txt','w')
out_fh.write('peptide_id'+'\t'+'pfam_id'+'\t'+'definition'+'\t'+'Evalues'+'\n')
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
 #   pfm_id=str(id_pfam[ids]).strip(']').strip('[')#you can split and join in future?
 #   pf_def=''#pfa_def means pfam[def]
 #   for p_ids in id_pfam[ids]:
 #       for defi in pfam_def[p_ids]:
 #           pf_def+=str(p_ids)+str([defi])+'; '
 #   pf_def=pf_def.strip().rstrip(';')
 #   pfm_id=pfm_id.strip(',')
 #   out_fh.write(str(ids)+'\t'+pfm_id+'\t'+pf_def+'\n')
#for ids in gid_pfam.keys():
#    print(ids)

