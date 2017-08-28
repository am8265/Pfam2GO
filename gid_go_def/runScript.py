from Pfam2GO import Pfam2GO


config = ""

pfam_dir = "/Users/ayanmalakar/BioinformaticsCookBook/Pfam_GO/gid_go_def/Pfam/Pfam-A.hmm"
protein_fasta_file = "/Users/ayanmalakar/BioinformaticsCookBook/first5seq.fa"
pfamtoGO_mapping_file = "/Users/ayanmalakar/BioinformaticsCookBook/Pfam_GO/gid_go_def/pfam2go"
go_runner = Pfam2GO(config)
output = go_runner.run_pfam_2_go(pfam_dir, protein_fasta_file, pfamtoGO_mapping_file)

#Testing..
print(output)