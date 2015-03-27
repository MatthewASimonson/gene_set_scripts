###################
# SNP Ratio Test:
###################

# Step 1: Generate SNP set files for each pathway

genes <- read.table("entrez.gene.map2",header=FALSE)
names(genes) <- c('CHR','START','STOP','GENE')
path <- read.table(list.files()[74],header=FALSE)
names(path) <- c('GENE','PATH','ANNOT')
gene.path <- merge(genes,path,by='GENE')
o.index <- order(gene.path$PATH)
gene.path <- gene.path[o.index,]

paths <- unique(gene.path$PATH)
path.gene.count <- vector() # count number of genes in each pathway
for(i in 1:length(paths)){# write out gene ranges for each pathway
index <- which(gene.path$PATH==paths[i])
path.gene.count[i] <- length(index)
path.set <- gene.path[index,c(2:4,1)]
path.name <- as.character(unique(gene.path$PATH[index]))
write.table(path.set,file=paste(path.name,".list",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
print(i)
}
#
for(i in 1:length(paths)){ # write out set files for each pathway
system(paste("plink --bfile hapmap-giant-common --make-set ",paths[i],".list --make-set-border 20 --write-set --out ",paths[i],"",sep=""))
}

# For pathways without top associations:
for (i in 1:length(paths)){
system(paste("plink --bfile hapmap-giant-notop --make-set ",paths[i],".list --make-set-border 20 --write-set --out notop.",paths[i],"",sep=""))
}

# Step 2: Create list with SNPs from each pathway and read in association data:

path.snps <- list()
path.snp.count <- vector()
for(i in 1:length(paths)){#
  system(paste("grep 'rs' ",paths[i],".set > ",paths[i],".snps",sep=""))# write list of SNPs from pathway to file
  snps <- read.table(paste("",paths[i],".snps",sep=""),header=FALSE)
  names(snps) <- c('SNP')
  path.snps[[i]] <- snps # assign pathway SNPs to list index
  path.snp.count[i] <- nrow(snps)
  print(i)
}
names(path.snps) <- paths

count.frame <- cbind.data.frame(paths,path.gene.count,path.snp.count)

# only look at pathways with between 20 and 200 genes:

sub.path.frame <- count.frame[which(path.gene.count>20 & path.gene.count<200),] 
sub.path.snps <- path.snps[as.numeric(row.names(sub.path.frame))]

snp.assoc <- read.table("GIANT.assoc",header=TRUE)
for(i in 1:length(sub.path.snps)){
  sub.path.snps[[i]] <- merge(sub.path.snps[[i]],snp.assoc,by="SNP")
  print(i)
}

save.image("SRT.Rdata")
