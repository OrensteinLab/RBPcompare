The scripts and functions


Data files:
eCLIP files:
	ENCODE_acession for example ENCFF933AGU
	Positive sequences file -	<ENCODE_acession>.bed.ext.hg19.fa
	Negative sequences file - 	<ENCODE_acession>.bed.control.ext.hg19.fa
	Structure files - 			<ENCODE_acession>.bed.ext.hg19.fa.combined_profile.txt.sho.txt
								Contain 5 rows for each sequence with the probabilities of each nucleotide to be paired/ hairpin- multi- inner- external-loop.
	
RNA BInd-n-Seq files:
	k-mer scores files -  		tab delimeted files containing a list of k-mers and their scores (as calculated by the authors) named {protein_name}_{score_type}_{k}.tsv (for example TIA1_Z_5.tsv).
	Structure files -			5 files for each experimet, each contain the probabilities for all nuclotide to be in a specific structure (calculated by RNAplfold)
	fastq files - 				containing all RNA sequences in the bound pool of the experiment

RNAcompete files:
	7-mer scores files -		Two files, one containing Z scores and the other E scores of all experimetns done by Ray 2013. each experiment is available in 3 columns, set A, set B and set AB.
	Sequence files - 			File containing all sequences with their messured binding intensity. 

Example raw data files are in the data directory

Supporting files:
	ProtANDConc.txt - 			Binary file containes a python dictionary with all RNA Bind-n-Seq relevant proteins and the concentrations in which the experiments were conducted in.
	eCLIP_BnS.txt -				Contain all paires of RBP and ENCODE identifier to eCLIP experiment in the overlap of the two data sets.
	BnS_CMPT.txt -				Contain the overlap of RNAcompete and RNA Bind-n-Seq. 
	CMPT_eCLIP.txt -			Contain the overlap of RNAcompete and eCLIP.
We provide all supporting files

kmers scores files:				We provide all 6mers scores and structure matricies as calculated by as. The scores are contained in binary files. 
								To load these file to python you can use the pickle package and use the following python commands:
										import pickle
										with open('file_name','rb') as fb:
											scores = pickle.load(fb)



installation:
there is no need to install the software but the dependencies specified below.

Dependencies:
python 3.6 and above
numpy
pandas
pickle
sklean
scipy


Code (python scripts):
CreateKmers-			Script for creating all possible k-mers.
						inputs: k, integer, the length of the RNA k-mer.
						
CountKmers-				This script will count the k-mers in the ".fastq" file.
						inputs (in that order): fastaqfile, k, Outfile.csv
						output: Outfile.csv - contain the counts for every k-mer
			
CountKmersWithProbs-	This script will count all k-mers and sum their probability matricies. 
						inputs (in that order): fastaqfile, k, E_profile_file,H_profile_file, I_profile_file, M_profile_file, Outfile.csv
						The output of this program will be a csv file, each row for one k-mer contains: counts and matix of structure probibilities:
						paired (alwayes will be zero and will be fixed in CreateScoresFiles.py), ______________.

eCLIP_BnS.py-			This script will assign a score to every RNA sequence in eCLIP experiment (both the positives and the control) and assign an AUC value based on the k-mers scores as computed by 
						RNA Bind-n-Seq developers.
						This script has no inputs as it conduct this process for all RBPs in eCLIP_BnS.txt file

BnS_CMPT.py- 			This script will use the RNA Bind-n-Seq score in the tsv files to predict binding from RNAcompete data.
						No inputs are required for the script since it will perfom that for all RBP-RNAcompete identifier paires in the BnS_CMPT.txt file/
						The output of this script is a csv file containing all combinations of concentration, frequency/ratio score, k and average/max aggregation function.
						the results will also be saved in a .txt format (as a binary file).

CMPT_CMPT.py-			This script will use the 7-mers scores from set A of RNAcompete experiment to predict binding of set B of the same RNAcompete experiment.
						It will be done for all RBPs in the overlap between RNAcompete and RNA BInd-n-Seq.
						The output is a csv file contain the usage of z/e scores and average/maximum aggregation function.

CMPT_eCLIP.py-			This script will assign a score to every RNA sequence in eCLIP experiment (both the positives and the control) and assign an AUC value based on the 7-mers scores as computed by 
						RNAcompete developers.
						This script has no inputs as it conduct this process for all RBPs in CMPT_eCLIP.txt file
						
CreateScoresFiles.py-	This script will fix the probabilities for pairedness and create a dictionary with the scores for each k-mer.
						inputs: protein - the name of the protein creating the score dictionary for.
								pup- use probability matrix as paried\unpaired or full 5 rows probabilities? 1- paired\unpaired 0- full matrix
								e-	use enrichment score (i.e ratio score)? 1-yes 0-no
								K-	length of k-mer, integer.
								m-	mean probabilities over nucleotide. 
						output: binary file containes the python dictionary of the results

GiveScores.py-			This script will asign a score for each sequence in the input file based on the scores file of the "CreateScoresFiles.py" script.
						Inputs:	protein- name of the protein to use it's scores.
								protein_to_score- the RNAcompete or ENCODE identifier of the protein. If the identifier starts with ENC it will be treated as an eCLIP identifier.
								pup- use probability matrix as paried\unpaired or full 5 rows probabilities? 1- paired\unpaired 0- full matrix
								e-	use enrichment score (i.e ratio score)? 1-yes 0-no
								K-	length of k-mer, integer.
								m-	mean probabilities over nucleotide. 						
						output- binary file that contain the AUC or Pearson correlation achived by the k-mers scores model. the name is ('identifier_e_pup_m.txt')