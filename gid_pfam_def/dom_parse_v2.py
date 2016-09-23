#!/usr/bin/env python
import re
import sys
import math
def gid_pfam(fh_1,eval=math.exp(-10)):
    
    fh=open(fh_1)
    gid_pfam=dict()
    pfam_def=dict()
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
            if evalue<=eval:#if there is a cutoff match....
                if pfam_def.get(pfam_id,'0')=='0':
                    pfam_def[pfam_id]=[def_pfam]#every unique pfam_id shoudl have a unique definition
                else:
                    pfam_def[pfam_id].append(def_pfam)
            
                if gid_pfam.get(gid,'0')=='0':
                    gid_pfam[gid]=[pfam_id]
                else:
                    gid_pfam[gid].append(pfam_id)
    
#print(gid_pfam)
    for ids in sorted(gid_pfam.keys()):
        gid_pfam[ids]=list(set(gid_pfam[ids]))
    for ids in pfam_def.keys():
        pfam_def[ids]=list(set(pfam_def[ids]))
    return (gid_pfam,pfam_def)
    fh.close()

gid_pfam=gid_pfam(sys.argv[1])
id_pfam=gid_pfam[0]
pfam_def=gid_pfam[1]
######testing here
print(id_pfam)
##########
out_fh=open('gid_pfam_def.txt','w')
out_fh.write('peptide_id'+'\t'+'pfam_id'+'\t'+'definition'+'\n')
for ids in id_pfam.keys():
    pfm_id=str(id_pfam[ids]).strip(']').strip('[')#you can split and join in future?
    pf_def=''#pfa_def means pfam[def]
    for p_ids in id_pfam[ids]:
        for defi in pfam_def[p_ids]:
            pf_def+=str(p_ids)+str([defi])+'; '
    pf_def=pf_def.strip().rstrip(';')
    pfm_id=pfm_id.strip(',')
    out_fh.write(str(ids)+'\t'+pfm_id+'\t'+pf_def+'\n')
#for ids in gid_pfam.keys():
#    print(ids)

