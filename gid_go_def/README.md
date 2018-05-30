# To follow a complete Workflow, use the 2 scripts listed below:


1) Pfam2GO.py

and runScript.py. 

Input Files

Pfam-A.hmm, protein fasta file, pfam2go file

Output files:

gid_to_go, pfam_to_go







# Additionally, Run the other Scripts like:

Run the Script parse_pfam2go_v2.py as: 

./parse_pfam2go_v2.py <pfam2go> <hmmscan_file>


Input Files:

pfam2go - maps Pfam Id to GO Id

hmmscan_file - output of hmmscan of query sequences in ‘domblout’ format

Output Files:

gid_go.txt - maps protein Id to GO Ids & GO Terms with Evalues




Parse_pfam2go_v2.py has the dependency:  
dom_parse.py script - parser script to help parse Hmmscan_file ( default evalue cutoff is 1e-10)

Run dom_parse.py as 
./dom_parse.py <hmmscan_file>




gid_pfam_def.txt - maps protein Id to PFAM Ids with PFAM definitions & PFAM Ids with e-values




E.g

./parse_pfam2go_v2.py pfam2go output1a.txt

./dom_parse.py output1a.txt

