gwas.dat <- read.table("MERGE.clean.FINAL.assoc.linear",header=TRUE)
genes <- read.table("clean.genelist.FINAL",header=TRUE)

gene.snp.ps <- vector()
gene.snp.rs <- vector()
for(i in 1:nrow(genes)){
  chr <- genes$CHR[i]
  start <- genes$START[i]
  stp <- genes$STOP[i]
  index <- which((gwas.dat$CHR==chr)&(gwas.dat$BP>start)&(gwas.dat$BP<stp))
  snp <- as.character(gwas.dat$SNP[index][which.min(gwas.dat$P[index])])
  if(length(snp)>0){
  gene.snp.rs[i] <- as.character(gwas.dat$SNP[index][which.min(gwas.dat$P[index])])
  gene.snp.ps[i] <- min(gwas.dat$P[index])
  }else{
  gene.snp.rs[i] <- NA
  gene.snp.ps[i] <- NA
  }  
  print(i)
}

gene.stats <- cbind.data.frame(genes[which(gene.snp.ps!='NA'),1],gene.snp.ps[which(gene.snp.ps!='NA')],gene.snp.rs[which(gene.snp.ps!='NA')])
names(gene.stats) <- c('GENE','P-VAL','SNP') # 12,687 genes have at least 1 SNP

# write out list of min SNP from each gene and calculate r^2 estimates using PLINK:
setwd("/STATGEN/home/simonsom/Obesity/Study_Merge/GWAS")
write.table(gene.stats[,3],file="min.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)
system("plink --bfile MERGE.clean.FINAL --extract min.snps --r2 --ld-window-r2 0 --out min.snp.r2")

# generate LD clumps:
assoc.total <- read.table("MERGE.clean.FINAL.assoc.linear",header=TRUE)
min.snps <- read.table("min.snps",header=FALSE)
names(min.snps) <- c('SNP')
ms.assoc <- merge(min.snps,assoc.total,by="SNP")
o.index <- order(ms.assoc$CHR,ms.assoc$BP)
ms.assoc <- ms.assoc[o.index,]
write.table(ms.assoc,file="ms.assoc.linear",quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')

#fix the linear file so it only contains min.snps
system("plink --bfile MERGE.clean.FINAL --extract min.snps --clump ms.assoc.linear --clump-p1 1 --clump-p2 1 --clump-r2 0.10 --clump-kb 250 --out min.snp.r2")
snp.clumps.pre <- read.table("min.snp.r2.clumped",header=TRUE)
snp.clumps <- snp.clumps.pre[!duplicated(snp.clumps.pre$SNP),]
snp.clumps <- snp.clumps[order(snp.clumps$CHR,snp.clumps$BP),]
gene.clumps <- merge(gene.stats,snp.clumps,by="SNP")[,c(1,2,3)]
clumped.stats <- gene.clumps[!duplicated(gene.clumps$SNP),]
write.table(clumped.stats,file="clumped.stats.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

# Count number of PCAs for each gene to correct for SNP multiple testing:
# Calculate PCA's for each gene:
################################
total.genes <- genes
names(total.genes) <- c('GENE','CHR','START','STOP')
rel.genes <- merge(clumped.stats,total.genes,by='GENE')

rel.gene.snps <- list() # fill with SNPs from each relevant gene
for(i in 1:nrow(rel.genes)){
  chr <- rel.genes$CHR[i]
  start <- rel.genes$START[i]
  stp <- rel.genes$STOP[i]
  index <- which((gwas.dat$CHR==chr)&(gwas.dat$BP>start)&(gwas.dat$BP<stp))
  snp <- as.character(gwas.dat$SNP[index])
  rel.gene.snps[[i]] <- snp
  print(i)
}

for(i in 1:nrow(rel.genes)){ # write list of snps for each relevant gene
write.table(rel.gene.snps[[i]],file=paste("",rel.genes$GENE[i],".snps",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
print(i)
}

rg.snp.count <- unlist(lapply(rel.gene.snps,length))
# Call multi-threaded packages
############################
CORES <- 14
require(foreach) # 
require(doMPI) # 
cl <- startMPIcluster(CORES) 
registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

loop_GENERATE_PED <- foreach(i=1:nrow(rel.genes)) %dopar% {
  gene <- as.character(rel.genes$GENE[i])
  system(paste("plink --bfile MERGE.clean.FINAL --noweb --extract ",gene,".snps --recode12 --out ",gene,"",sep="")) # recode data as non-binary and remove subjects without SNPs
 } # end loop_GENERATE_PED
#

 keep_PCA <- function(fit,R2){ # return number of PCA's that explain 'R2' percent of variance
    var.PC <- as.numeric(fit$sdev)^2
    total.var <- sum(var.PC)
    greater.index <- which(cumsum(var.PC)>R2*total.var)
    keep.count <- min(greater.index)
    return(keep.count)
  }

correct.p <- vector()
for(i in 1:nrow(gene.PCA.stats.pre)){
correct.p[i] <- p.adjust(gene.PCA.stats.pre[i,3], method = "bonferroni", n = gene.PCA.stats.pre[i,4])
print(i)
}
gene.stats.uc <- gene.stats # assign uncorrected gene.stats

gene.stats.c <- cbind.data.frame(gene.PCA.stats.pre$GENE,gene.PCA.stats.pre$SNP,correct.p)
names(gene.stats.c) <- c('GENE','SNP','CORRECTED-PVAL')
gene.stats.c <- gene.stats.c[!duplicated(gene.stats.c$SNP),]
gene.clumps <- merge(gene.stats.c,snp.clumps,by="SNP")
gene.clumps <- gene.clumps[order(gene.clumps$CHR,gene.clumps$BP),]
plot(-1*log10(gene.stats[,2]))
gene.final <- gene.clumps[,c(2,3)]
write.table(gene.stats.c[,c(1,3)],file="c.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')


# Examine P-value for SNPs in top pathways:

sig.paths <- read.table("Sig_Pathways.txt",header=TRUE,sep="\t",fill=NA)

sig.list <- as.list(sig.paths)

merge_gene <- function(sig.list,gene.stats){
path.data <- as.data.frame(sig.list)
names(path.data) <- c('GENE')
path.stats <- merge(path.data,gene.stats,by='GENE')
return(path.stats)
}

total.path.stats <- lapply(sig.list,merge_gene,gene.stats)

path.lengths <- lapply(total.path.stats,nrow)

extract_p <- function(total.path.stats){
return(sort(-1*log10(total.path.stats[,2])))
}

path.pvals <- lapply(total.path.stats,extract_p)
names(path.pvals) <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17')

boxplot(path.pvals, ylab="-log10(p)",xlab="Pathway",main="Top Gene Associations Grouped By Pathway")
save.image("snp.analysis")

# Examine Gene Overlap Between Pathways:

overlap.matrix <- matrix(NA,nrow=17,ncol=17)

for(i in 1:ncol(sig.paths)){
  for(j in 1:ncol(sig.paths)){
  overlap.matrix[i,j] <- length(intersect(sig.paths[,i],sig.paths[,j]))-1
  } # end loop j
} # end loop i

# Probability of observed overlap:
hyper.matrix <- matrix(NA,nrow=17,ncol=17)

for(i in 1:ncol(sig.paths)){
  for(j in 1:ncol(sig.paths)){
    d.white <- overlap.matrix[i,j] # number of genes in pathway 1 and pathway 2; overlap
    v.white <- length(unique(c(sig.paths[,i],sig.paths[,j])))-1  # total number of genes in pathway 1 or pathway 2
    d.black <- v.white-d.white # number of genes in pathway 1 or pathway 2 that are not in pathway 1 and pathway 2
    v.black <- 6804-v.white 
  hyper.matrix[i,j] <- phyper(d.white,v.white,v.black,d.white+d.black,lower.tail=FALSE)
  } # end loop j
} # end loop i


# Correspondence Analysis:
percent.matrix <- matrix(NA,nrow=17,ncol=17)

for(i in 1:ncol(sig.paths)){
  for(j in 1:ncol(sig.paths)){
  percent.matrix[i,j] <- (length(intersect(sig.paths[,i],sig.paths[,j]))-1)/(length(unique(c(sig.paths[,i],sig.paths[,j])))-1)
  } # end loop j
} # end loop i

sig.in <- c(2,3,4,8,10,11,15) # index with replicated pathways

raw.dist <- as.data.frame(1-percent.matrix[sig.in,sig.in])
names(raw.dist) <- c('KEGG Long Term Depression','Reactome Axon Guidance','KEGG Calcium Signaling','KEGG Axon Guidance','Reactome TACS','KEGG O-Glycan Biosynthesis','KEGG GRNH Signaling','KEGG Focal Adhension','Reactome Neuro','Reactome Regulation of Insulin Secretion GLP1','KEGG ARVC','KEGG Complement and Coagnulation C','Reactome Formation and M','Reactome Centrosome M','KEGG Endocytosis','Reactome Formation of T','KEGG Base Excision Repair')[sig.in]
dist.matrix <- as.dist(raw.dist)
fit <- hclust(dist.matrix,method="ward")
plot(fit, main="Clustering by Genes Shared",xlab="Pathways That Replicate Across Methods")
groups <- cutree(fit, k=3) # cut tree into 3 clusters

# How often do genes appear?

non.axon.genes <- as.character(unlist(sig.list[c(1,3,5:17)])) # in all other pathways besides axon guidance
axon.genes <- unique(as.character(unlist(sig.list[c(2,4)])))

path.genes <- c(axon.genes,non.axon.genes)

gene.keep.index <- which(table(path.genes)<100) # remove empty space
gene.tab <- table(path.genes)[gene.keep.index]

# examine if overlap count is related to association strength

path.gca <- cbind.data.frame(names(gene.tab),gene.tab) # path names and overlap count
names(path.gca) <- c('GENE','COUNT')

stat.gca <- merge(path.gca,gene.stats,by="GENE")

overlap.mod <- lm(-log10(stat.gca[,3])~stat.gca[,2])
boxplot(-log10(stat.gca[,3])~stat.gca[,2], xlab="Number Of Times Gene Appears In Top Pathways", ylab="-log10(p)",main="Shared Genes Influence on Obesity")
