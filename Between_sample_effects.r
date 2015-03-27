# Examine between sample bias that may be occuring:
###################################################

# First load data:
setwd("/home/simonsom/Obesity/Study_Merge/GenePath")

dataset <- 'MERGE.clean.FINAL'
#load("Cross_Validation.Rdata")
load("resid.models.Rdata") # load total phenotype with covariates for all subjects
total <- total[(!is.na(total$IID) & !is.na(total$BMI) & !is.na(total$SEX) & !is.na(total$AGE) & !is.na(total$BATCH) & !is.na(total$C1)),] # remove any subjects with missing data

# Read in list of all genes that have PCAs:
####################################################################################################
#system("ls *.min_PCA | sed 's/\\(.*\\)......../\\1/' > PCA.genes") # write list of all genes that have .min_PCA format to file
PCA.genes <- read.table("PCA.genes",header=FALSE)
names(PCA.genes) <- c('NAME')
gene.count <- nrow(PCA.genes)

ARIC.index <- which(total$set=='ARIC')
CARDIA.index <- which(total$set=='CARDIA')
MESA.index <- which(total$set=='MESA')
  
  ARIC.mod <- lm(formula = as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + 
    as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + 
    as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + 
    as.numeric(C9) + as.numeric(C10) + as.numeric(C11) + as.numeric(C12) + 
    as.numeric(C13) + as.numeric(C14) + as.numeric(C15) + as.numeric(C16) + 
    as.numeric(C17) + as.numeric(C18) + as.numeric(C19) + as.numeric(C20) + as.factor(BATCH), data = total[ARIC.index,])

  CARDIA.mod <- lm(formula = as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + 
    as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + 
    as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + 
    as.numeric(C9) + as.numeric(C10) + as.numeric(C11) + as.numeric(C12) + 
    as.numeric(C13) + as.numeric(C14) + as.numeric(C15) + as.numeric(C16) + 
    as.numeric(C17) + as.numeric(C18) + as.numeric(C19) + as.numeric(C20) + as.factor(BATCH), data = total[CARDIA.index,])

  MESA.mod <- lm(formula = as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + 
    as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + 
    as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + 
    as.numeric(C9) + as.numeric(C10) + as.numeric(C11) + as.numeric(C12) + 
    as.numeric(C13) + as.numeric(C14) + as.numeric(C15) + as.numeric(C16) + 
    as.numeric(C17) + as.numeric(C18) + as.numeric(C19) + as.numeric(C20) + as.factor(BATCH), data = total[MESA.index,])

  total.mod <- lm(formula = as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + 
    as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + 
    as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + 
    as.numeric(C9) + as.numeric(C10) + as.numeric(C11) + as.numeric(C12) + 
    as.numeric(C13) + as.numeric(C14) + as.numeric(C15) + as.numeric(C16) + 
    as.numeric(C17) + as.numeric(C18) + as.numeric(C19) + as.numeric(C20) + as.factor(set) + as.factor(BATCH), data = total)

  ARIC.phe<- resid(ARIC.mod)
  MESA.phe<- resid(MESA.mod)
  CARDIA.phe<- resid(CARDIA.mod)

  total.phe <- resid(total.mod)

  ARIC.set <- cbind.data.frame(rep(0,length(ARIC.phe)),total[ARIC.index,1],ARIC.phe)
  names(ARIC.set) <- c('FID','IID','PHE')

  MESA.set <- cbind.data.frame(rep(0,length(MESA.phe)),total[MESA.index,1],MESA.phe)
  names(MESA.set) <- c('FID','IID','PHE')

  CARDIA.set <- cbind.data.frame(rep(0,length(CARDIA.phe)),total[CARDIA.index,1],CARDIA.phe)
  names(CARDIA.set) <- c('FID','IID','PHE')

ARIC.stat <- as.matrix(unlist(lapply(gene.PCAs,lin.mod,ARIC.set)))
MESA.stat <- as.matrix(unlist(lapply(gene.PCAs,lin.mod,MESA.set)))
CARDIA.stat <- as.matrix(unlist(lapply(gene.PCAs,lin.mod,CARDIA.set)))


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

  stats <- ARIC.stat[,1] # assign z-value from each gene to array

  ES.ARIC=getES(set.idx=set.idx,gene.stats=stats) # returns deviation from what is expected by random chance given normal distribution (Enrichment Scores)

  NES.ARIC=(ES-mES)/sdES # normalized observed enrichment scores for observed data; ONE FOR EACH PATHWAY

  stats <- MESA.stat[,1] # assign z-value from each gene to array

  ES.MESA=getES(set.idx=set.idx,gene.stats=stats) # returns deviation from what is expected by random chance given normal distribution (Enrichment Scores)

  NES.MESA=(ES.MESA-mES)/sdES # normalized observed enrichme

  stats <- CARDIA.stat[,1] # assign z-value from each gene to array

  ES.CARDIA=getES(set.idx=set.idx,gene.stats=stats) # returns deviation from what is expected by random chance given normal distribution (Enrichment Scores)

  NES.CARDIA=(ES.CARDIA-mES)/sdES # normalized observed enrichme


plot(tot.NES$NES~tot.NES$Path)
abline(aric.mod,col="green")
abline(cardia.mod,col="blue")
abline(mesa.mod,col="red")
