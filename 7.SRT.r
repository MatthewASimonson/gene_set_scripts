############################################################
# FOLLOWING WERE MODIFIED BY MATTHEW A. SIMONSON 2011-2012 #
############################################################
##
## the following functions are for running ODushlaine et al. (2009) method.

###################
# SNP Ratio Test: #
###################

# INPUT:
# 1.) name of dataset
# 2.) number of cores
# 3.) number of permutations
# 4.) name of file of list of unique autosomal genes with header (NAME, START, STOP, CHR) ;
# 5.) threshold of interest for SRT
# 6.) +- intergenic region surrounding gene to include as genic in KB
# NOTE: A .RDATA FILE WITH PATHWAY GENE LOCATIONS MUST BE PROVIDED: generate_gene_set.r will generate this file
# NOTE: A PHENOTYPE FILE MUST BE PROVIDED: format of 2 columns (Ind ID, Pheno) no header, 'datasetname.pheno' file name
# NOTE: A COVARIATE FILE MUST BE PROVIDED: format of 24 columns (Ind ID, Sex, Age, Batch, PCA1, PCA2, ... , PCA20) no header, 'datasetname.covar' file name

# OUTPUT:
#1.)

SRT <- function(dataset,CORES,PERMS,clean.genelist,threshold,KB){ # START SRT UMBRELLA FUNCTION
  
  generate_snp_path <- function(dataset,clean.genelist,KB){ # This function should be run first to generate the list of SNPs within pathways
  total.snps <- read.table(paste('',dataset,'.bim',sep="")) # read in all SNPs
  names(total.snps) <- c('CHR','SNP','CM','BP','A1','A2')
  
  gene.locations <- read.table(clean.genelist,header=TRUE) # read in all genes

  snp.info.data <- as.data.frame(cbind(as.character(total.snps$SNP),total.snps$CHR,total.snps$BP))
  gene.info.data <- as.data.frame(cbind(as.character(gene.locations$NAME),gene.locations$CHR,gene.locations$START,gene.locations$STOP))
  #
 Snp2Gene <- function(snp.info, gene.info, dist=KB*1000){ # Map SNPs to Genes
  # INPUTS:
  # snp.info -> like a map file; 3 cols, 1) snp name 2) CHR 3) base position
  # gene.info -> like a map file for genes; 4 cols, 1) gene name 2) CHR 3) Start bp 4) Stop bp
  # dist -> number of bp up and down stream of gene boundaries to consider genic
  ################################################################################################
    snp.info <- as.matrix(snp.info)
    gene.info <- as.matrix(gene.info)
    Chr.nam <- unique(as.vector(gene.info[,2])) # create list of chromosomes included in list of genes
    gene.snplist <- list()
    gene.snplist[[nrow(gene.info)+1]] <- 0 # fill last index of list with zero
    snps.pos <- as.numeric(snp.info[,3]) # create object to store base positions of each SNP
    snps.name <- snp.info[,1] # store SNP ID's
    
    start <- as.integer(as.vector(gene.info[,3])) # store all gene start position
    end <- as.integer(as.vector(gene.info[,4])) # store all gene stop positions
    gene.pos <- cbind(start,end) # create matrix of all start and end positions
    gene.rang <- cbind(gene.pos[,1]-dist, gene.pos[,2]+dist) # add dist to range of gene boundaries
       
    for (i in Chr.nam){
       idx.gene <- which(gene.info[,2]==i)
       idx.snp <- which(snp.info[,2]==i)
       snps.pos.i <- snps.pos[idx.snp] # get SNP bp's for current chromosome
       snps.name.i <- snps.name[idx.snp] # get list of SNPs for current chromosome

       for (j in idx.gene){ # loop through indeces in gene file for current chromosome
           gene.rang.j <- gene.rang[j,] # assign the range of each genic region
           snps.j <- snps.name.i[which(snps.pos.i>=gene.rang.j[1]&snps.pos.i<=gene.rang.j[2])] # index of SNP bp's in current chromosome with within genic range
           if (length(snps.j)>0){
             gene.snplist[[j]] <- snps.j # if some snps exist within a gene, then assign them to list of each gene
           } # END if
         } # END LOOP j
       } # END LOOP i
       gene.snplist <- gene.snplist[1:nrow(gene.info)] # remove zero index of gene.snplist
       names(gene.snplist) <- gene.info[,1] # give proper names to list of genes
       return(gene.snplist)  
     } # END FUNCTION Snp2Gene

    Gene2Set <- function(gene, set){ # 
    set.idx <- lapply(set, function(x){
        idx <- NULL    
        for (i in x) idx <- c(idx, which(gene==i)) # loop through each position in pathway and find index with name of gene from list
        return(idx)
        })
    return(set.idx)
    } # END FUNCTION 'Gene2Set'

     transform <- function(gene.set){
      values <- as.matrix(gene.set[[1]])[,1] # Transform gene-pathway list into correct format
      return(values)
    } # END FUNCTION 'transform'
  
  Snp2Path <- function(set.idx,gene.snps){ # insert SNPs from each gene into their respective pathways
    for(j in 1:length(set.idx)){
       values <- set.idx[[j]]
       fill <- vector()
       for (i in 1:length(values)){
         fill <- c(gene.snps[[values[i]]],fill)
       }# END LOOP i     
       set.idx[[j]] <- fill
    } # END loop j    
    snp.paths <- set.idx # apply function 'snp_fill' to every position in set.idx
     # Convert list into data.frame format:
     snp.path.matrix <- matrix(ncol=2)
     for(i in 1:length(snp.paths)){
       path.snps <- snp.paths[[i]]
       path.names <- rep(names(snp.paths[i]),length(path.snps))
       col.bind <- cbind(path.names,path.snps)
       snp.path.matrix <- rbind(col.bind,snp.path.matrix)
     } # END LOOP i
     snp.path.frame <- as.data.frame(snp.path.matrix[1:(nrow(snp.path.matrix)-1),])
    return(snp.path.frame) 
  } # END FUNCTION 'Snp2Path

    gene.snps <- Snp2Gene(snp.info.data,gene.info.data,KB) # takes about 30 seconds to run for whole genome
    load(paste('',dataset,'.gene.set.Rdata',sep="")) # read in gene set information created using 'generate_gene_set()' function
    gene.set <- gene.sets.final
    gene.sets <- lapply(gene.set,transform)    
    set.idx <- Gene2Set(gene.info.data[,1], gene.sets) # indeces in list of genes in pathway sets
    snp.paths <- Snp2Path(set.idx,gene.snps)

 write.table(snp.paths,file=paste("",dataset,".SRT_paths.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t') # Also write out with correct name
} # END FUNCTION generate_snp_path
  
  perm.assoc <- function(dataset,CORES,PERMS){ # permute phenotype information and generate plink .assoc files for each permutation and non-permuted

  ############################  
  # Read in covariate data:  #
  ############################

  fam.data <- read.table(paste(dataset,'.fam',sep=""),header=FALSE)
  names(fam.data) <- c('FID','IID')
  covar.data <- read.table(paste(dataset,'.covar',sep=""),header=FALSE) # read in covariate data
  names(covar.data) <- c('IID')
  covar.data <- covar.data[!duplicated(covar.data$IID),] # remove any duplicates
  fam.ind <- fam.data[,c(1,2)]
  m.covar <- merge(fam.ind,covar.data,by='IID')
  covar.assoc <- m.covar[,c(2,1,5,3,4,6:ncol(m.covar))] 
  names(covar.assoc) <- c('FID','IID','BATCH','SEX','AGE','PCA1','PCA2','PCA3','PCA4','PCA5','PCA6','PCA7','PCA8','PCA9','PCA10','PCA11','PCA12','PCA13','PCA14','PCA15','PCA16','PCA17','PCA18','PCA19','PCA20')

  ############################  
  # Read in phenotype data:  #
  ############################
    
  phe.data <- read.table(paste(dataset,'.phe',sep=""),header=FALSE)
  names(phe.data) <- c('IID')
  m.fam.data <- merge(covar.assoc[,1:2],phe.data,by='IID')
  assoc.phe <- m.fam.data[,c(2,1,3)]
  names(assoc.phe) <- c('FID','IID','BMI')

  ###########################################  
  # Merge data and residualize covariates:  #
  ###########################################

  phe.covar <- na.omit(merge(assoc.phe,covar.assoc,by='IID')) # merge and remove NA's
  phe.covar <- phe.covar[!duplicated(phe.covar$IID),] # remove any duplicates

  covar.model <- lm(BMI ~ BATCH + SEX + AGE + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + PCA19 + PCA20, data=phe.covar)

  resid.BMI <- resid(covar.model) # get part of obesity not predicted by covariates

  resid.phe <- as.data.frame(phe.covar$FID.x)
  resid.phe$IID <- phe.covar$IID
  resid.phe$BMI <- resid.BMI
  names(resid.phe) <- c('FID','IID','BMI')

  ###########################################  
  # Write out residualized phenotype data:  #
  ###########################################

  write.table(resid.phe,file=paste(dataset,'.assoc_pheno',sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE)# write out unpermuted

  for (i in 1:PERMS){ # generate permuted phenotypes
    assoc.phe.perm <- resid.phe
    assoc.phe.perm[,3] <- sample(resid.phe[,3],replace=FALSE) # permute phenotype data (sampling without replacement)
    write.table(assoc.phe.perm,file=paste(dataset,'.assoc_pheno',i,'',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE) # write out phenotype file for initial association to be used by Plink
  } # END LOOP i

  ##########################################################################################
  # Modify Plink format files to contain only SNPs within specified pathways of interest:  #
  ##########################################################################################
  snp.paths <- read.table(paste('',dataset,'.SRT_paths.txt',sep=""),header=FALSE)
  write.table(snp.paths$V2,file="snps.in.paths",row.names=FALSE,col.names=FALSE,quote=FALSE) # write out list of SNPs that are within pathways being examined
  
  system(paste('plink --bfile ',dataset,' --extract snps.in.paths --make-bed --out ',dataset,'.GENE',sep=""))# create new Plink format files with .GENE in name
  
  #################################
  # Call multi-threaded packages: #
  #################################

  require(foreach) # 
  require(doMPI) # 
  cl <- startMPIcluster(CORES) 
  registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

  #####################################################################################
  # perform association on unpermuted and permuted data in using multi-threaded loop: #
  #####################################################################################

  system(paste('nohup plink --bfile ',dataset,'.GENE --pheno ',dataset,'.assoc_pheno --linear --out ',dataset,' &',sep=""))  # for unpermuted 

  # Permutations:
  loop_assoc <- foreach(j=1:PERMS) %dopar% {
    system(paste('plink --bfile ',dataset,'.GENE --pheno ',dataset,'.assoc_pheno',j,' --linear --out ',dataset,'.perm',j,' ',sep=""))
  } # END LOOP j

  #
} # END FUNCTION perm.assoc


parse.sig.runSRT <- function(dataset,threshold,PERMS){
###########################################
# Prepare association files and read in:  #
###########################################
  
 system("ls *.linear > linear.temp")
  system("ls *.linear | sed 's/\\(.*\\)......./\\1/' > assoc.temp") # remove .linear extension
  temp.linear <- read.table("linear.temp",header=FALSE)
  temp.assoc <- read.table("assoc.temp",header=FALSE)
  for (i in 1:nrow(temp.linear)){
   system(paste('cp ',temp.linear[i,1],' ',temp.assoc[i,],'',sep="")) # copy all .linear files to .assoc format
 }
  system("rm linear.temp; rm assoc.temp") # remove temp files

  system("perl parse_assoc_files.pl") # convert all .assoc files to .forSRT format (SNPid,p-value)

# Read in association files

observed.snp.p <- read.table(paste(dataset,'.assoc.forSRT',sep=""),header=TRUE) 
perm.data <- as.data.frame(matrix(NA,nrow=nrow(observed.snp.p),ncol=PERMS+2))

# Store observed associations and SNP name in first 2 columns:
perm.data[,1] <- observed.snp.p$SNP
perm.data[,2] <- observed.snp.p$P
names(perm.data) <- c('SNP','ObservedP')

# Read in permuted data:
for(h in 1:PERMS){
print(paste("reading permutation",h,"",sep=""))
perm.data[,h+2] <- read.table(paste(dataset,'.perm',h,'.assoc.forSRT',sep=""),header=TRUE,colClasses=c('character','numeric'),nrow=nrow(observed.snp.p))$P
}

#############################
# Read in pathway with SNP: #
#############################
 
snp.path <- read.table(paste(dataset,'.SRT_paths.txt',sep=""),header=FALSE)
names(snp.path) <- c('PATH','SNP')
path.perms <- merge(snp.path,perm.data,by='SNP') # merge each SNPs p-vals with their pathwway label
paths <- unique(path.perms$PATH)

########################################################################################################
# Find SNP ratio for each pathway and empirical significance for that pathway BEFORE multiple testing: #
########################################################################################################
 
path.data <- as.data.frame(matrix(0,nrow=length(paths),ncol=PERMS+2))
names(path.data) <- c('PATH','Observed.RATIO')

 SNP.count <- vector() # store size of each pathway

 for(i in 1:length(paths)){
    temp.path <- path.perms[which(path.perms$PATH==paths[i]),3:ncol(path.perms)]
    SNP.count[i] <- nrow(temp.path)
 }

 count.frame <- cbind.data.frame(paths,SNP.count) # store size of each pathway with label
 
   SR_row <- function(temp.path,cutoff){
     sig.count <- length(which(temp.path<=cutoff))
     total.count <- length(temp.path)
     SR <- sig.count/total.count
     return(SR)
   }

 for (j in 1:length(paths)){ # loop through all pathways
     print(paste('Examining pathway: ',j,'',sep=''))
     cutoff <- threshold
     temp.path <- path.perms[which(path.perms$PATH==paths[j]),3:ncol(path.perms)]
     path.data[j,2:ncol(path.data)] <- apply(temp.path,2,SR_row,cutoff)
   } # end loop j
 #
   
path.data$PATH <- paths
 
save(path.data,file=paste('',dataset,'.path.data.Rdata',sep="")) # save object to current directory 

# Calculate p-values from SNP Ratios:
 
 get_pval <- function(ob,null){
  pct.a<- mean(null >= ob)
  M <- length(null)
  B1 <- pct.a*M
  pval <- (B1+1)/(M+1)
  return(pval)
  }

p.path.data <- path.data
names(p.path.data) <- c('PATH','Observed.P')
for(i in 1:nrow(p.path.data)){
  SRs <- jitter(as.vector(path.data[i,3:ncol(p.path.data)],'numeric'),2)
  SR.ob <- path.data[i,2]
  p.path.data[i,3:ncol(p.path.data)] <- sapply(SRs,get_pval,SRs)
  p.path.data[i,2] <- get_pval(SR.ob,SRs)
}

uc.pvals <- p.path.data[,2]
write.table(uc.pvals,file=paste(dataset,'.uc.SRT',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE) 
                  
# Find corrected p-values:

 #########################################################
 # FDR based on empirical null distribution of p-values: #
 #########################################################
 
 mp.dist <- as.matrix(p.path.data[,3:ncol(p.path.data)])# null p-values
 null.dist <- jitter(mp.dist,2) # add a little noise to make data less discrete

 q.vals <- vector()
 for(i in 1:length(uc.pvals)){
   nulls <- mean(null.dist<=uc.pvals[i])*length(null.dist)
   obs <- mean(uc.pvals<=uc.pvals[i])*length(null.dist)
   q.vals[i] <- (nulls+1)/(obs+1) # predicted number of false discoveries per number of observed
 }


 
 # Store p-values:

 path.pvals <- as.data.frame(paths)
 names(path.pvals) <- "PATHWAY"

 path.pvals$'P-VAL' <- round(uc.pvals,6)
 path.pvals$'Q-VAL' <- round(q.vals,6)
 o.index <- order(path.pvals[,2])
 final.pvals <- path.pvals[o.index,]
 
write.table(final.pvals,file=paste(dataset,'.SRT',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
#
# 
} # END FUNCTION parse.applySRT


##########################
# CALL NESTED FUNCTIONS: #
##########################
  
 generate_snp_path(dataset,clean.genelist,KB) # write out pathway information in correct format
 perm.assoc(dataset,CORES,PERMS) # Run GWAS and permute data
 parse.sig.runSRT(dataset,threshold,PERMS) # Execute SRT on data and assess pathway significance (which is adjusted for multiple testing); NOTE: THIS STEP MAY TAKE SEVERAL HOURS DUE TO PERMUTATION TESTING

  
} # END FUNCTION SRT

################################
# RUN SRT UMBRELLA FUNCTION:   #
################################

dataset <- 'CARDIA.clean.FINAL'
CORES <- 5
PERMS <- 1000
clean.genelist <- "clean.genelist.FINAL"
threshold <- 0.01
KB <- 20

SRT(dataset,CORES,PERMS,clean.genelist,threshold,KB) # START SRT UMBRELLA FUNCTION


resid.phe <- read.table("MERGE.clean.FINAL.assoc_pheno",header=FALSE)
