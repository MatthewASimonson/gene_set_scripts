# Function Inputs:
# 1. name of dataset #
# 2. r^2 >= cutoff for variance explained by PCA's from original data EXPRESSED AS DECIMAL (85% = .85)
#
# Function Outputs:
# 1. GRM for each gene in dataset
# 2. .eigenval and .eigenvec file for each gene
# 3. .min_PCA file for each gene in dataset

# NOTE: limit this to 3 cores per node on RC, could use up to 3Gb RAM per cores

generate_gene_PCA <- function(dataset,R2,CORES,clean.list){ # dataset must be in quotes, BE SURE R2 IS EXPRESSED AS DECIMAL, example: .85

# Switch to directory containing binary plink files:
####################################################
curdir <- getwd()
setwd(paste('',curdir,'/',dataset,'_plink_bin',sep=""))
  
# Read in list of all genes that have binary plink files (if they contain any SNPs on current array):
####################################################################################################
   no_snps <- function(set.data){ # RETURN LIST OF GENES THAT HAVE NO SNPS
      no.snp.genes <- vector()
      end.index <- which(set.data$V1=="END")
      pre.index <- (end.index-1) # last SNP (if any) in gene
      has.snps <- set.data[pre.index,] %in% grep("rs", set.data[pre.index,], value = TRUE)
      if(length(unique(has.snps))==1 & unique(has.snps)==FALSE){
        has.snps <- set.data[pre.index,] %in% grep("SNP_A", set.data[pre.index,], value = TRUE) # If SNPs are affy instead of rs
     }
      no.snp.genes <- set.data[pre.index,][has.snps==FALSE]
      drop.genes <- as.data.frame(no.snp.genes)
      return(drop.genes)
    } # END FUNCTION 'no_snps'

set.data <- read.table(paste('',dataset,'.set',sep=""),header=FALSE) # read in .set file
gene.list <- read.table(clean.genelist,header=TRUE) # read in generated gene list
drop.genes <- no_snps(set.data)
drop.gene.index <- match(drop.genes$no.snp.genes,gene.list$NAME)    

keep.genes <- gene.list[c(-drop.gene.index),] # new list of genes that have SNPs and thus plink files
bin.genes <- as.data.frame(keep.genes$NAME)
names(bin.genes) <- 'NAME'
gene.count <- nrow(bin.genes)

# Calculate PCA's for each gene:
################################

 keep_PCA <- function(fit,R2){ # return number of PCA's that explain 'R2' percent of variance
    var.PC <- as.numeric(fit$sdev)^2
    total.var <- sum(var.PC)
    greater.index <- which(cumsum(var.PC)>R2*total.var)
    keep.count <- min(greater.index)
    return(keep.count)
  }

# Call multi-threaded packages
############################

require(foreach) # 
require(doMPI) # 
cl <- startMPIcluster(CORES) 
registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

loop_GENERATE_PED <- foreach(i=1:nrow(bin.genes)) %dopar% {
  gene <- as.character(bin.genes$NAME[i])
  system(paste("plink --bfile ",gene," --noweb --recode12 --out ",gene,"",sep="")) # recode data as non-binary and remove subjects without SNPs
 } # end loop_GENERATE_PED

loop_MIN_PCA <- foreach(i=1:nrow(bin.genes)) %dopar% {
  gene <- as.character(bin.genes$NAME[i])
  data <- read.table(paste("",gene,".ped",sep=""),header=FALSE)
  # Calculate sum of alleles for each loci: 4,3,2=AA,Aa,aa
  geno.mat.c1 <- data[,7:ncol(data)]
  geno.mat.c2 <- cbind.data.frame(data[,8:ncol(data)],rep(0,nrow(data)))
  geno.pre <- geno.mat.c1 + geno.mat.c2
  col.index <- seq(1,ncol(geno.pre),2) # only keep odd columns
  geno.matrix <- geno.pre[,col.index]
  fit <- princomp(geno.matrix,cor=TRUE)
  keep.count <- keep_PCA(fit,R2)
  total.PCA <- fit$scores
  min.PCA.pre <- total.PCA[,1:keep.count] # only keep PCA's that explain specified percent of variance for given gene
  min.PCA <- cbind.data.frame(data[,1:2],min.PCA.pre) # merge FID and IID column
  write.table(min.PCA,file=paste("",gene,".min_PCA",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
} 

#
} # end function


dataset <- 'MESA.clean.FINAL'
R2 <- .85
CORES <- 6
clean.genelist <- 'clean.genelist.FINAL'

generate_gene_PCA(dataset,R2,CORES,clean.list)



