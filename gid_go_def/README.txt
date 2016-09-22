Run the Script like:

./parse_pfam2go_v2.py <annotation_file from pfam2go> <output of hmmscan of query sequences in -domblout format>

eg. :

./parse_pfam2go_v2.py pfam2go output1a.txt

Result:

gid_go.txt file would be generated mapping peptide_ids to GO terms and GO id information and GO_id,Evalue information

Please Note the dependencies: 1)dom_parse.py script 
