#Here we use an example code to show the drug repositioning steps as follows:
#1.Construct a similarity matrix;
#2.Optimize the similarity matrix by screening the optimal drug dose/time duration; 
#3.Optimize the similarity matrix by screening the optimal shRNA perturbagen;
#4.Rank the drugs by absolute correlation coefficients (ACC) and extract the top 10 drugs with highest ACCs as the final drugs.

#We provided three folders:
#Test: the example data, A549 and HCC515 "SH" and "CP" transcriptomics profiles (truncated gctx files, level 5), saved in the /Test/sh and /Test/cp directories,respectively. Folders Output1/2/3 were used to save the output results from the these steps.

#meta_data: siginfo_beta_yuan.Rdata is same as the "siginfo_beta.txt" which was downloaded from the LINCS data portal (https://clue.io/data/CMap2020#LINCS2020) ;  drug_info.txt is a processed file based on the "compoundinfo_beta.txt", in which we extracted the drugs' pert_id and drug commonly used names.

#Program: two main programs used for the prediction.