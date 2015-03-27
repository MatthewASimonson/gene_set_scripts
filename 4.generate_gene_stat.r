# Function Inputs:
# 1. name of dataset
# 2. number of cores
# 3. number of permutations
# NOTE: PLINK FORMAT FILES FOR THE DATASET MUST BE IN WORKING DIRECTORY
# NOTE: A PHENOTYPE FILE MUST BE PROVIDED: format of 2 columns (Ind ID, Pheno) no header, 'datasetname.pheno' file name
# NOTE: A COVARIATE FILE MUST BE PROVIDED: format of 24 columns (Ind ID, Sex, Age, Batch, PCA1, PCA2, ... , PCA20) no header, 'datasetname.covar' file name
#
# Function Outputs:
# 1. single row temporary file for each gene ('.stat_temp' extension)
# 2. file with one column of genes and one column of z-values ('.gene_stat' extension)


generate_gene_stats <- function(dataset,CORES,perms){ # inputs must be character (in quotes)

# Read in list of all genes that have PCAs:
####################################################################################################
system("ls *.min_PCA | sed 's/\\(.*\\)......../\\1/' > PCA.genes") # write list of all genes that have .min_PCA format to file
PCA.genes <- read.table("PCA.genes",header=FALSE)
names(PCA.genes) <- c('NAME')
gene.count <- nrow(PCA.genes)

  ############################  
  # Read in covariate data:  #
  ############################

#  fam.data <- read.table(paste(dataset,'.fam',sep=""),header=FALSE)
#  names(fam.data) <- c('FID','IID')
#  covar.data <- read.table(paste(dataset,'.covar',sep=""),header=TRUE) # read in covariate data
#  names(covar.data) <- c('IID')
#  covar.data <- covar.data[!duplicated(covar.data$IID),] # remove any duplicates
#  names(covar.data) <- c('IID')
#  fam.ind <- fam.data[,c(1,2)]
#  m.covar <- merge(fam.ind,covar.data,by='IID')
#  covar.assoc <- m.covar[,c(2,1,5,3,4,6:ncol(m.covar))] 
#  names(covar.assoc) <- c('FID','IID','BATCH','SEX','AGE','PCA1','PCA2','PCA3','PCA4','PCA5','PCA6','PCA7','PCA8','PCA9','PCA10','PCA11','PCA12','PCA13','PCA14','PCA15','PCA16','PCA17','PCA18','PCA19','PCA20')

  ############################  
  # Read in phenotype data:  #
  ############################
    
 # phe.data <- read.table(paste(dataset,'.pheno',sep=""),header=FALSE)
 # names(phe.data) <- c('IID')
 # m.fam.data <- merge(covar.assoc[,1:2],phe.data,by='IID')
 # assoc.phe <- m.fam.data[,c(2,1,3)]
  

  ###########################################  
  # Merge data and residualize covariates:  #
  ###########################################

#  phe.covar <- na.omit(merge(assoc.phe,covar.assoc,by='IID')) # merge and remove NA's
#  phe.covar <- phe.covar[!duplicated(phe.covar$IID),] # remove any duplicates

#  covar.model <- lm(PHE ~ BATCH + SEX + AGE + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + PCA19 + PCA20, data=phe.covar)

#  resid.PHE <- resid(covar.model) # get part of phenotype not predicted by covariates

#  resid.phe <- as.data.frame(phe.covar$FID.x)
#  resid.phe$IID <- phe.covar$IID
#  resid.phe$PHE <- resid.PHE
  resid.phe <- read.table(paste(dataset,'.phe',sep=""),header=FALSE)
  names(resid.phe) <- c('FID','IID','PHE')

# Read in eigenSNPs for each gene and place dataframes in list
##############################################################

read_gene <- function(gene){ # Function that reads gene listed in data.frame of gene names (1 per row)
     gene.PCAs.temp <- read.table(paste("",gene,".min_PCA",sep=""),header=FALSE) # read in eigenSNPs for gene
     gene.PCAs.temp <- gene.PCAs.temp[,2:ncol(gene.PCAs.temp)]
     eigen.names <- paste(rep('EIG',(ncol(gene.PCAs.temp)-1)),1:(ncol(gene.PCAs.temp)-1),sep="")
     names(gene.PCAs.temp) <- c('IID',eigen.names)
     gene.PCAs.temp <- as.matrix(gene.PCAs.temp)
     return(gene.PCAs.temp)
   }

gene.PCAs <- apply(PCA.genes, 1, read_gene) # apply read_gene function to PCA.genes list of genes data.frame
names(gene.PCAs) <- PCA.genes$NAME

lin.mod <- function(gene.PCA,resid.phe){ # Function that return z-value for each gene
  gene.covars <-  merge(gene.PCA,resid.phe,by='IID')#
  eigen.SNP.text <- c(paste(rep('EIG',(ncol(gene.PCA)-1)),1:(ncol(gene.PCA)-1),sep=""))
  eigen.SNP.variable <- paste(eigen.SNP.text,collapse=" + ")
  model.A <- eval(parse(text=(paste("lm(PHE ~ ",eigen.SNP.variable,", data=gene.covars)",sep="")))) # evaluate effects of each eigenSNP after controlling for covariates
  pval <- as.numeric(pf(summary(model.A)$fstatistic[1],summary(model.A)$fstatistic[2],summary(model.A)$fstatistic[3],lower.tail=FALSE))
  zval <- -1*qnorm(pval/2) # store z-value
return(zval)
}

# Generate un-permuted phenotpe gene statistics:

gene.stat <- as.data.frame(unlist(lapply(gene.PCAs,lin.mod,resid.phe)))

gene.stat <- cbind.data.frame(PCA.genes,gene.stat)
names(gene.stat) <- c('GENE','Z-VAL')
write.table(gene.stat,file=paste('',dataset,'.gene_stat',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)

# Start permutation loop:
#########################

# Call multi-threaded packages

require(foreach) # 
require(doMPI) # 
cl <- startMPIcluster(CORES) 
registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

# START MPI/FOREACH LOOP
loop <- foreach(i=1:perms) %dopar% {
# first permute phenotype:
  perm.phe <- as.matrix(resid.phe)
  perm.phe[,3] <- sample(perm.phe[,3],replace=FALSE)
  
gene.perm.stat <- as.data.frame(unlist(lapply(gene.PCAs,lin.mod,perm.phe))) 
names(gene.perm.stat) <- c('p.z')
# write out each permutation:
write.table(gene.perm.stat,file=paste('',dataset,'.perm_stat',i,'',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
} # END MPI/FOREACH LOOP
#
# Read in temporary files for each gene and combine into single file then write out
################################################################################### 

gm2 <- matrix(0,nrow=gene.count,ncol=(perms)) # permutation matrix to be filled
fill <- as.data.frame(gm2) # data frame to be filled with permutations
for (i in 1:perms){
 fill[,i] <-  read.table(paste('',dataset,'.perm_stat',i,'',sep=""),header=TRUE)
 print(i)
} # END FUNCTION
fill <- cbind.data.frame(PCA.genes,fill) # add label column
write.table(fill,paste('',dataset,'.gene_perms',sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
#
} # END FUNCTION generate_gene_stats

dataset <- 'CARDIA.clean.FINAL'
clean.genelist <- 'clean.genelist.FINAL'
perms <- 1000
CORES <- 6
