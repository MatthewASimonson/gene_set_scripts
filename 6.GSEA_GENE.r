######################################################################################################################
# FOLLOWING WERE MODIFIED (AND LARGELY WRITTEN) BY MATTHEW A. SIMONSON 2011-2012
######################################################################################################################
## the following two functions are the ones for running Wang et al. (2007) method.


getES <- function(set.idx, gene.stats, p=1){ # This function is called in the next function; returns deviation from what is expected by random chance given normal distribution for each pathway
    ES=rep(0,length(set.idx)) 
    rk <- rank(-gene.stats,ties.method="first") # return order index based on ranked z-values (first of a tie is chosen as ordered first); greatest to least
    N=length(gene.stats)  ## total number of genes

    for (i in 1:length(set.idx)){
      path.idx=set.idx[[i]]  ## the gene indice for this i-th pathway
      Nh=length(path.idx)  ## number of genes in this pathway
      oo <- sort(rk[path.idx]) # sort ranks of genes for path i
  
      ES.all= -(1:N)/(N-Nh) # 1 through total number of genes, divided by total genes - genes in pathway
      statj=gene.stats[path.idx]
      statj=-sort(-statj)
      Nr=sum(abs(statj)^p) #
      for (j in 1:(Nh-1)){ # loop through number of genes in pathway
         jj=sum(abs(statj[1:j])^p)
         ES.all[oo[j]:(oo[j+1]-1)] = ES.all[oo[j]:(oo[j+1]-1)]+jj/Nr+j/(N-Nh)
      }   
      ES.all[N]=0
      ES[i]=max(ES.all)
    }
    return(ES)
  }


getp.gsea <- function(set.idx, stats, stats0, path.names){ # This function is called at the end of gseaSNP;  inputs 1.) genes in each pathway list, 2.) true stats for each gene, 3.) permuted stats for each gene

  ES=getES(set.idx=set.idx,gene.stats=stats) # returns deviation from what is expected by random chance given normal distribution (Enrichment Scores)

  ES0=NULL
  for (i in 1:ncol(stats0)){ 
    s0=stats0[,i] 
     ES0=cbind(ES0,getES(set.idx,gene.stats=s0)) # get null Enrichment Scores for each permutation of data
  }

  mES=apply(ES0,1,mean) # find the mean Enrichment Scores of null
  sdES=apply(ES0,1,sd) # find standard deviation of null Enrichment Scores
  NES=(ES-mES)/sdES # normalized observed enrichment scores for observed data; ONE FOR EACH PATHWAY

  NES0 <- (ES0-mES)/sdES # normalized null

# Calculate uncorrected p-values: 

 get_pval <- function(ob,null){
  pct.a<- mean(null >= ob)
  M <- length(null)
  B1 <- pct.a*M
  pval <- (B1+1)/(M+1)
  return(pval)
  }
  
 uc.pvals <- vector() # store uncorrected p-values

 uc.pvals <- sapply(NES,get_pval,NES0)

 pvals[which(pvals[,1]=='KEGG_FOCAL_ADHESION'),]

 #########################################################
 # FDR based on empirical null distribution of p-values: #
 #########################################################

#null.dist <- 2*pnorm(-1*abs(NES0))

# q.vals <- vector()
# for(i in 1:length(uc.pvals)){
#   nulls <- mean(null.dist<=uc.pvals[i])*length(null.dist)
#   obs <- mean(uc.pvals<=uc.pvals[i])*length(null.dist)
#   q.vals[i] <- (nulls+1)/(obs+1) # predicted number of false discoveries per number of observed
# }

 # Store p-values:
 
  pval <- cbind.data.frame(path.names,round(uc.pvals,6))
  pvals <- cbind.data.frame(path.names,ES)
  names(pvals) <- c('PATH','ES')
  return(pvals)
}


# INPUTS:
# 1.) dataset = name of data set in character vector
# 2.) gene.info = name of gene list file

# OUTPUTS:
# 1.) array of p-values for each gene based on 
gseaGene <- function(dataset,gene.info.file,perms){

    # read in data:
    z.data <- read.table(paste('',dataset,'.snp_stat',sep=""),header=TRUE)# read in z-values for each gene, names of indeces should be gene names
    z.permutes <- read.table(paste('',dataset,'.gene_perms',sep=""))## read in matrix with permutation values
    gene.names <- as.data.frame(z.data[,1]) # names of all genes with z scores
    names(gene.names) <- c('NAME')
    gene.info.pre <- read.table("clean.genelist.FINAL",header=TRUE) # read in list of genes
    # order gene list correctly:
    gene.info <- merge(gene.names,gene.info.pre,by='NAME')
    
    load(paste('',dataset,'.gene.set.Rdata',sep="")) # read in gene set information
   #
    
    Gene2Set <- function(gene, set){
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
    
    gene.sets <- lapply(gene.sets.final,transform)    
    set.idx <- Gene2Set(gene.info[,1], gene.sets) # indeces in list of genes of genes in pathway sets 
    #
    stats <- z.data[,2] # assign z-value from each gene to array
    stats0 <- z.permutes[,2:ncol(z.permutes)]# assign z-value from each gene for each perm to column of dat.frame
    path.names <- names(gene.sets)
    pvals <- getp.gsea(set.idx, stats, stats0, path.names) # now use above functions and permutation test results to return pvals for each pathway; input= 1.) genes in each pathway list, 2.) true stats for each gene, 3.) permuted stats for each gene, 4.) names of pathways in same order as stats
    write.table(pvals,file=paste('',dataset,'.GSEA',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
}

