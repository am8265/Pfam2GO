#!/usr/bin/env python
import sys
import re
from dom_parse import gid_pfam  # gid_pfam function is imported from dom_parse module 
def gid_go(fh1, fh2):
    '''pfam2go(<pfam2go file>,<hmmscan_output>) and return value is a dictionary of dictionaries ----gid_go_eval={gid:{(go,go_term):[evals],.....}}'''
    fh1 = open(fh1)  # pfam2GO file
    gid_pfam_eval = gid_pfam(fh2)  # gid_pfam from hmmscan--output1.txt for pfam domains-output files-  gid -pfam, also have pfam to definition!
# Looking into pfam2GO file-----------    
    pfam2go = dict()
    for line in fh1:
        if re.search('^Pfam', line) != None:
            line = line.rstrip('\n').strip()
            parse_line = line.split(':')[1].split(' ')
            pfam = parse_line[0]
            # print(pfam)
            parse_line = ':'.join(line.split(':')[2:])
            go = parse_line
            if pfam2go.get(pfam, 0) == 0:  # if pfam id doesn't exist then create a new one
                # pfam2go[pfam]=go
                 pfam2go[pfam] = [go]
            # print(pfam2go)
            else:  # if the pfam id already exist as a key then append GO Terms to previous list
# pfam2go[pfam]=pfam2go[pfam]+','+go
                pfam2go[pfam].append(go)
                pfam2go[pfam] = list(set(pfam2go[pfam]))
    print(pfam2go)
    fh1.close()
# Making the gid_go_eval--------------
    gid_go_eval = dict()
# Looking into Output1a.txt file for pepide_id to Pfam_id informtion----------------
    for gid in sorted(gid_pfam_eval.keys()):
        print(gid)
        gid_go_eval[gid] = {}  # everytime we encounter a UNIQUE gene id we create a new dictionary from it ...a test of gid_go_eval[gid]={} if for that gid no GO ANNOTATION exists for ALL pfam_ids found on it!
        for pfams,evals in gid_pfam_eval[gid].items():
            pfam_id = pfams[0]  # get me only the pfam_id!
            if pfam2go.get(pfam_id, 'NA') != 'NA':  # pfam2GO exist for this pfam_id!
                gos=pfam2go.get(pfam_id)  # note GO refers to GO and GO_term
                for go in gos:
                # parsing go to get GO and GO_term
                    if go.find(';') != -1:
                        go = go.split(';')
                        go_def = go[0].strip()
                        go = go[1].strip()                                                                                                                       
                    else:
                        go_def = 'NA'  # if in certain small cases GO_definition isn't there for a GO accesion 
                    if gid_go_eval[gid].get((go,go_def), 'NA') == 'NA':
                        gid_go_eval[gid][(go,go_def)] = evals
                    else:
                        gid_go_eval[gid][(go,go_def)].extend(evals)
            else:
                pass
    return gid_go_eval
gid_go_eval = gid_go(sys.argv[1], sys.argv[2])  # sys.argv[1] should be the pfam2go file, sys.argv[2] is the output1.txt file from hmmscan 
# print(gid_go)
# gid_pfam=gid_go_pfam[1]
fh = open('gid_go.txt', 'w')
# out_pfam=open('gid_pfam_def.txt','w')
fh.write('peptide_id'+'\t'+'GO_id'+'\t'+'GO_Def'+'\t'+'Evalues'+'\n')  # change this line for changing the columns 
for id in sorted(gid_go_eval.keys()):
    cat_go = ''  # variable for catenating GO ids                                                                                                                               
    cat_go_def = ''  # variable for catenating GO def and GO term                                                                                                               
    cat_go_eval = ''
    if gid_go_eval[id] != {}:  # i.e there is atleast 1 GO annotation for this peptide_id!
        for go_inf,evals in gid_go_eval[id].items():
            print(go_inf)
            print(evals)
            go = go_inf[0]
            go_def = go_inf[1]
            min_eval = min(evals)
            cat_go += go + ','
            cat_go_def += go + '[' + go_def + ']; '
            cat_go_eval += go + '[' + str(min_eval) + ']' + '; '
        # Processin the FINAL LINE for a gid!
        cat_go = cat_go.rstrip(',')
        cat_go_def = cat_go_def.rstrip().rstrip(';')
        cat_go_eval = cat_go_eval.rstrip().rstrip(';')
        fh.write(id + '\t' + cat_go + '\t' + cat_go_def + '\t' + cat_go_eval + '\n')
    
    else:
         fh.write(id + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\n')  # for this peptide_id NO GO ANNOTATION EXIST!
fh.close()    

