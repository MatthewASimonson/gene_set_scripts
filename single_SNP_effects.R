setwd("/STATGEN/home/simonsom/Obesity/Study_Merge/Sig_Paths")

# Examine the associations of SNPs within pathways that are significant in the total sample
sp <- read.table("Sig_Path_Overlap.txt",header=TRUE,sep='\t',fill=TRUE)
sp.mat <- as.matrix(sp)

table(sp.mat) # 

overlap.mat <- matrix(NA,nrow=5,ncol=5)

RAG.KCS <- intersect(sp.mat[,1],sp.mat[,2]) # 1,91,6804-91,100; phyper=0.3880284

RAG.KAG <- intersect(sp.mat[,1],sp.mat[,3]) # 31,91,6804-91,130; phyper=p<.0001
RAG.RRI <- intersect(sp.mat[,1],sp.mat[,4]) # 0
RAG.KAR <- intersect(sp.mat[,1],sp.mat[,5]) # 2,91,6804-91,76; phyper=0.08106437
KCS.KAG <- intersect(sp.mat[,2],sp.mat[,3]) # 1,100,6804-100,130; phyper=0.5736364
KCS.RRI <- intersect(sp.mat[,2],sp.mat[,4]) # 5
KCS.KAR <- intersect(sp.mat[,2],sp.mat[,5]) # 1
KAG.KAR <- intersect(sp.mat[,3],sp.mat[,4]) # 2
RRI.KAR <- intersect(sp.mat[,4],sp.mat[,5]) # 5
