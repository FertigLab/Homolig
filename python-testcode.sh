#homolig_format_checker.py confirms format of the input file. It will report to you any uninterpretable V gene calls and the number of duplicated or otherwise invalid amino acid sequences. 

#homolig_wrapper.py calculates similarity scores between sequences using the similarity matrix specified in --metric. 
# pairwise mode will calculate all similarity scores among members resulting in nxn output matrix in .h5ad format. 
# axb mode will cross two input files resulting in axb output matrix. See homolig_wrapper.py -h

#homolig_write_umap.py clusters the resultant similarity matrix from homolig_wrapper.py. 
# It will also report silhouette scores in output file.

#homolig_write_umap.py generates UMAP coordinates based on the input similiarity matrix. 	

#h5ad_to_csv.py converts .h5ad output format of homolig_wrapper.py into a .csv file for those interested in manipulating the raw data. 	

#Three series of examples are provided in this code. 
#(i) Paired chain data from anti-viral TCRs in VDJDb used in manuscript. 
#(ii) Example TCRB data from one of the HC samples used in manuscript. 
#(iii) Example IGH data from one of the HC samples used in manuscript.  

####Process paired-chain viral TCR data


python3 ./homolig/homolig_format_checker.py\ 
		-i ./test-data/TCR-paired-test_viral.csv\
		--seq tcr \
		--chains paired \
		--metric ensemble \
		--mode pairwise \
		--save_germline False\
		--verbose True \
		--output ./vdjdb_paired_ens.h5ad


python3 ./homolig/homolig_wrapper.py \
		-i ./test-data/TCR-paired-test_viral.csv\
		--seq tcr \
		--chains paired \
		--metric ensemble \
		--mode pairwise \
		--save_germline False\
		--verbose True \
		--output ./vdjdb_paired_ens.h5ad
		
 		
python3 ./homolig/clusterHomolig.v4.py \
		-i ./vdjdb_paired_ens.h5ad \
		--num_clusters 80

python3 ./homolig/homolig_write_umap.py \
	-i ./vdjdb_paired_ens.h5ad\
	-o ./vdjdb_paired_ens_umap.csv	
	
python3	./homolig/h5ad_to_csv.py \
	-i ./vdjdb_paired_ens.h5ad
	
	
####Process single-chain TCRB data

python3 ./homolig/homolig_format_checker.py\
		-i ./test-data/TCRB-test-Subject_12.csv\
		--seq tcr \
		--chains beta \
		--metric ensemble \
		--mode pairwise \
		--save_germline False\
		--verbose True \
		--output ./TCRB_ens.h5ad

python3 ./homolig/homolig_wrapper.py \
		-i ./test-data/TCRB-test-Subject_12.csv\
		--seq tcr \
		--chains beta \
		--metric ensemble \
		--mode pairwise \
		--save_germline False\
		--verbose True \
		--output ./TCRB_ens.h5ad
		
		
python3 ./homolig/clusterHomolig.v4.py \
		-i ./TCRB_ens.h5ad \
		--num_clusters 80
		
python3 ./homolig/homolig_write_umap.py \
	-i ./TCRB_ens.h5ad \
	-o ./TCRB_ens_umap.csv	
	
python3	./homolig/h5ad_to_csv.py \
	-i ./TCRB_ens.h5ad 
	
	
####Process single-chain IGH data
python3 ./homolig/homolig_format_checker.py\
		-i ./test-data/IGH-test_D1-M_0_BRR_D1-M-001.csv\
		--seq bcr \
		--chains heavy \
		--metric pam250 \
		--mode pairwise \
		--save_germline False\
		--verbose True \
		--output ./IGH_pam.h5ad

python3 ./homolig/homolig_wrapper.py \
		-i ./test-data/IGH-test_D1-M_0_BRR_D1-M-001.csv\
		--seq bcr \
		--chains heavy \
		--metric pam250 \
		--mode pairwise \
		--save_germline False\
		--verbose True \
		--output ./IGH_pam.h5ad
		
		
python3 ./homolig/clusterHomolig.v4.py \
		-i ./IGH_pam.h5ad \
		--num_clusters 80
		
python3 ./homolig/homolig_write_umap.py \
	-i ./IGH_pam.h5ad \
	-o ./IGH_pam_umap.csv	
	
python3	./homolig/h5ad_to_csv.py \
	-i ./IGH_pam.h5ad 
