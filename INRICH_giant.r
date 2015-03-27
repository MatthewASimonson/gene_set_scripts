# Generate INRICH input files:
################################
path.count <- 1452 # number of pathways in file
# create target gene set file with entrez gene ids:

gene.set.list <- list()

for(i in 0:(path.count-1)){ # read in each pathway
gene.set.list[[i+1]] <- read.table("c2.cp.v3.1.entrez.gmt",header=FALSE,skip=i,nrow=1)
print(i)
}

gene.set.frame.list <- list()
for(i in 1:path.count){ # convert each row into data frame
path.frame <- cbind.data.frame(t(gene.set.list[[i]][3:length(gene.set.list[[i]])]),rep(unlist(gene.set.list[[i]][2]),length(gene.set.list[[i]])-2),rep(unlist(gene.set.list[[i]][1]),length(gene.set.list[[i]])-2))
names(path.frame) <- c('GENE_ID','GENE_Set','ANNOTATION')
gene.set.frame.list[[i]] <- path.frame
print(i)
}

gene.set.frame <- gene.set.frame.list[[1]] # merge all pathways into single data frame
for(i in 2:path.count){
gene.set.frame <- rbind(gene.set.frame,gene.set.frame.list[[i]])
print(i)
}
write.table(gene.set.frame,file="msigbd3.0.set",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

