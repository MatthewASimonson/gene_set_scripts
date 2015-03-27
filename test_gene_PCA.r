
dataset <- 'ARIC2'
R2 <- 85
CORES <- 1

genes <- read.table("temp.list",header=FALSE)
bin.genes <- as.data.frame(genes$V4)
names(bin.genes) <- c('NAME')
gene.count <- nrow(bin.genes)

# Read in number of SNPs in each gene and append column to gene data frame:
###########################################################################
system('rm genes.snp_count') # remove previous version of this file
for (i in 1:gene.count){
system(paste("grep 'SNPs selected from' ",bin.genes$NAME[i],".log|cut -d ' ' -f 1 >> genes.snp_count",sep=""))
}
SNP.count <- read.table("genes.snp_count",header=FALSE) # read in gene
names(SNP.count) <- 'SNPs'
bin.genes <- cbind(bin.genes,SNP.count)
#

# Call multi-threaded packages:
###############################

require(foreach) # 
require(doMPI) # 
cl <- startMPIcluster(CORES) 
registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

# Make GRM and perform eigen-decomposition for each gene using multi-threaded loop:
###################################################################################

loop_GRM <- foreach(i=1:gene.count) %dopar%

for i in (1:gene.count){
  gene <- as.character(bin.genes$NAME[i])    
  system(paste("gcta --bfile ",gene," --make-grm --save-ram --out ",gene,"",sep="")) # generate GRM (Run time depends on number of SNPs in gene)
  system(paste("gcta --grm ",gene," --pca ",bin.genes$SNPs[i]," --out ",gene,"",sep="")) # eigen-decomposition
}

# START HERE
# Read in .eigenval files and calculate number of eigenvectors that explain >= R2% of SNP variance for each gene:
#################################################################################################################

loop_R2 <- foreach(i=1:gene.count) %dopar% {
  gene <- as.character(bin.genes$NAME[i])  
  gene.eigval <- abs(as.matrix(read.table(paste("",gene,".eigenval",sep=""),header=FALSE))) # read in eigenvals for gene

  total.var <- sum(as.matrix(gene.eigval)) # total variance
  num.PCAs <- 0 # create counter for number of PCA required to predict R2% of SNP variance
  accum.R2 <- 0
  PRE <- 0
  j <- 1
  while (PRE<R2){
    j <- j+1
    num.PCAs <- j
    accum.R2 <- sum(as.matrix(gene.eigval[1:j]))
    PRE <- accum.R2/total.var
    print(j)
  } # end while loop that calculates how many PCA's required

system(paste("cut -d ' ' -f 1-",num.PCAs+2," ",gene,".eigenvec > ",gene,".min_PCA",sep=""))  
} # end loop_R2
