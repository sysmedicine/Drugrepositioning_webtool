#Main code
#1.Construct similarity matrices for A549 and HCC515
source("../Program/DrugRepositioning_Functions.R")
cell_list<-as.matrix(c("A549","HCC515"))
path_input<-"../Test/"
path_output<-paste0(path_input,"Output_1/")
result<-DR_cal_corr(path_input,path_output,cell_list)
##########################################################
#2.Optimize the similarity matrix by screening the optimal drug dose/time duration; 
#It is needed to filter the positive and negative correlation coefficients,respectively. Here we only show how to filter positive correlation coefficients.
rm(list=ls())
source("../Program/DrugRepositioning_Functions.R")
cell_list<-as.matrix(c("A549","HCC515"))
dir="pos"#"pos" or "neg"
load("../meta_data/siginfo_beta_yuan.Rdata")
path_input="../Test/Output_1"
path_output="../Test/Output_2"
for (i in 1:length(cell_list)){
  result<-DR_filter_drug(path_input,path_output,sig_info,dir,cell_list[i])
}
##########################################################
#3.Optimize the similarity matrix by screening the optimal shRNA perturbagen;
rm(list=ls())
source("../Program/DrugRepositioning_Functions.R")
cell_list<-as.matrix(c("A549","HCC515"))
dir="pos"#"pos" or "neg"
load("../meta_data/siginfo_beta_yuan.Rdata")
path_input="../Test/Output_2"
path_output="../Test/Output_3"
for (i in 1:length(cell_list)){
  result=DR_filter_other_trt(path_input,path_output,sig_info,dir,cell_list[i])
}
##########################################################
#4.Rank the drugs by the decreasin order of absolute correlation coefficients (ACC) and extract the top 10 drugs with highest ACCs as the final drugs.
#Example: extract top 10 by selecing A549 cell line
rm(list=ls())
int_target="MCM6"#query gene
int_cell="A549"
#effect_type<-"inhibitor"#inhibitor or activator
#pert_type<-"sh"#perturbation type of data
abs_rcoef=0.4#absolute correlation coefficient no less than 0.4 
top_num_drugs=10#extract the top 10 drugs

drug_info<-as.data.frame(read.csv(file="../meta_data/drug_info.txt",header=T,sep="\t"))
path_raw<-"../Test/Output_3/"
setwd(path_raw)
df<-readRDS(paste0(int_cell,"_Rmatrix_FilterBoth.RDS"))
df<-as.matrix(df[,int_target])
colnames(df)<-int_target
loc<-match(rownames(df),drug_info$pert_id)
ext_drug<-cbind(as.data.frame(drug_info[loc,]),as.data.frame(df))

ext_drug<-ext_drug[which(ext_drug[,3]>=abs_rcoef),]#ACC>=0.4
rank<-order(ext_drug[,3],decreasing = T)
ext_drug<-ext_drug[rank,]
ext_drug<-ext_drug[1:top_num_drugs,]
ext_drug
