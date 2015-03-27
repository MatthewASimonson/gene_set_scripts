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


# Start permutation loop where 1000 random 50% splits of the data are created and evaluated:
#########################

perms <- 100
# START FOR LOOP

for(i in 1:perms){ # iterate through perms
# Designate 50% of sample to training set and 50% to test set:
  
  training.index <- sample(1:nrow(total),(nrow(total)/2),replace=FALSE)
  training.mod <- lm(formula = as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + 
    as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + 
    as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + 
    as.numeric(C9) + as.numeric(C10) + as.numeric(C11) + as.numeric(C12) + 
    as.numeric(C13) + as.numeric(C14) + as.numeric(C15) + as.numeric(C16) + 
    as.numeric(C17) + as.numeric(C18) + as.numeric(C19) + as.numeric(C20) + 
    as.factor(set) + as.factor(BATCH), data = total[training.index,])
  resid.phe.train <- resid(training.mod)
  
  test.mod <- lm(formula = as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + 
    as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + 
    as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + 
    as.numeric(C9) + as.numeric(C10) + as.numeric(C11) + as.numeric(C12) + 
    as.numeric(C13) + as.numeric(C14) + as.numeric(C15) + as.numeric(C16) + 
    as.numeric(C17) + as.numeric(C18) + as.numeric(C19) + as.numeric(C20) + 
    as.factor(set) + as.factor(BATCH), data = total[c(-training.index),])
  resid.phe.test <- resid(test.mod)
  
  training.set <- cbind.data.frame(rep(0,length(resid.phe.train)),total[training.index,1],resid.phe.train)
  names(training.set) <- c('FID','IID','PHE')
  
  test.set <- cbind.data.frame(rep(0,length(resid.phe.test)),total[c(-training.index),1],resid.phe.test)
  names(test.set) <- c('FID','IID','PHE')
  
training.stat <- as.matrix(unlist(lapply(gene.PCAs,lin.mod,training.set)))
test.stat <- as.matrix(unlist(lapply(gene.PCAs,lin.mod,test.set)))

names(training.stat) <- c('p.z')
names(test.stat) <- c('p.z')
# write out each permutation:
write.table(training.stat,file=paste('',dataset,'.cv_train_stat',i,'',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
write.table(test.stat,file=paste('',dataset,'.cv_test_stat',i,'',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
  #
} # END LOOP i
#
# Read in temporary files for training and test set permutations
################################################################

test.stats <- as.data.frame(matrix(0,nrow=7161,ncol=(perms))) # matrix to be filled
training.stats <- as.data.frame(matrix(0,nrow=7161,ncol=(perms))) # matrix to be filled

for (i in 1:perms){
 test.stats[,i] <-  read.table(paste('',dataset,'.cv_test_stat',i,'',sep=""),header=TRUE)
 training.stats[,i] <-  read.table(paste('',dataset,'.cv_train_stat',i,'',sep=""),header=TRUE)
 print(i)
} 

# Run GSEA on training set and test set data
############################################

    z.permutes <- read.table(paste('',dataset,'.gene_perms',sep=""))## read in matrix with permutation values
    gene.names <- as.data.frame(PCA.genes$NAME) # names of all genes 
    names(gene.names) <- c('NAME')
    gene.info.pre <- read.table('clean.genelist.FINAL',header=TRUE) # read in list of genes
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
    stats0 <- z.permutes[,2:ncol(z.permutes)]# assign z-value from each gene for each perm to column of dat.frame
    path.names <- names(gene.sets)

################################################
# Get p-values for training sets and test sets #
################################################

getES <- function(gene.stats,set.idx, p=1){ # This function is called in the next function; returns deviation from what is expected by random chance given normal distribution for each pathway
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

  ES.test <- apply(test.stats,2,getES,set.idx) # Enrichment Scores for test sets
  ES.train <- apply(training.stats,2,getES,set.idx) # Enrichment Scores for training sets

  ES0=NULL
  for (i in 1:ncol(stats0)){ 
    s0=stats0[,i] 
     ES0=cbind(ES0,getES(set.idx,gene.stats=s0)) # get null Enrichment Scores for null data
  }

  mES=apply(ES0,1,mean) # find the mean Enrichment Scores of null
  sdES=apply(ES0,1,sd) # find standard deviation of null Enrichment Scores

  NES.test=(ES.test-mES)/sdES # normalized observed enrichment scores for observed data; ONE FOR EACH PATHWAY
  NES.train=(ES.train-mES)/sdES

  NES0 <- (ES0-mES)/sdES # normalized null

# Calculate uncorrected p-values: 

 get_pval <- function(ob,null){
  pct.a<- mean(null >= ob)
  M <- length(null)
  B1 <- pct.a*M
  pval <- (B1+1)/(M+1)
  return(pval)
  }
  
 uc.pvals.train <- matrix(0,nrow=512,ncol=perms) # store uncorrected p-values
 uc.pvals.test <- matrix(0,nrow=512,ncol=perms)

cross_val_p <- function(data){  
return(sapply(data,get_pval,NES0))
}

uc.pvals.train <- apply(NES.train,2,cross_val_p)
uc.pvals.test <- apply(NES.test,2,cross_val_p)

save.image(file="Cross_Validation.Rdata")
# Examine the correlation across cross validation samples of p-values from pathways
#########################################################################################

below_per <- function(train,test,cutoff,path.dat){ # returns % of times pathway has p-value below cutoff in both training set AND test set, also p-value for percentage based on empirical null
 train.below <- train 
 train.below[which(train<=cutoff)] <- 1
 train.below[which(train>cutoff)] <- 0 
 
 test.below <- test
 test.below[which(test<=cutoff)] <- 1
 test.below[which(test>cutoff)] <- 0

 total.below <- test.below + train.below
 total.below[which(total.below<2)] <- 0
 total.below[which(total.below==2)] <- 2
 per.2 <- function(x){return(mean(x==2))}
 percent.overlap <- apply(total.below,1,per.2)
 paths.below <- cbind.data.frame(path.dat,percent.overlap)
 names(paths.below) <- c('PATHWAY','PERCENT')
 return(paths.below) # what % of the time is data equal to or less than specified cutoff
}
cutoff <- .05
dat <- below_per(uc.pvals.train,uc.pvals.test,cutoff,path.names)

pval.data <- read.table("MERGE.clean.FINAL.GSEA",header=TRUE)
names(pval.data) <- c('PATHWAY','P.VAL','Q.VAL')
sig.paths <- pval.data[which(pval.data[,2] <= cutoff ),1]
path.names <- names(gene.sets)
sig.path.index <- match(sig.paths,path.names)

null.paths <- pval.data[which(pval.data[,2] > cutoff ),1]
path.names <- names(gene.sets)
null.path.index <- match(null.paths,path.names)
top.paths <- dat[sig.path.index,]
bottom.paths <- dat[null.path.index,]
t.test(top.paths[,2],bottom.paths[,2]) # compare match rates in those above and below cutoff
#

# Create Plot of pathway p-values:

total.paths.pre <- merge(rbind(top.paths,bottom.paths),pval.data,by='PATHWAY')
tot.pval.index <- order(total.paths.pre[,3])
total.paths <- cbind.data.frame(total.paths.pre[tot.pval.index,],1:512)
plot(total.paths[,2],xlab="Pathway Ranked By Total Sample P-Value",ylab="% Training and Test Both Below Cutoff",main="GSEA Pathway Cross-Validation",pch=19,bg='grey')
p <- .05/512
null.05 <- mean(uc.pvals.train<=.05)*mean(uc.pvals.test<=.05) # probability any path cross validation replicates by chance
bon.p <- qexp(p,1/null.05,lower.tail=FALSE) # .05 cutoff before bonferroni
abline(h=bon.p,col="red",lty=5,lwd=2)

legend(300,.2, # places a legend at the appropriate place
c("Positive","Type I","Type II","Negative","p = 0.05"), # puts text in the legend
pch=c(19,19,19,19,NA),       
lty=c(NA,NA,NA,NA,5), # gives the legend appropriate symbols (lines)
lwd=c(NA,NA,NA,NA,2),       
col=c("blue","green3","orange","black","red"))

points(total.paths[((total.paths[,3]<=.05) & (total.paths[,2]>=bon.p)),5],total.paths[((total.paths[,3]<=.05) & (total.paths[,2]>=bon.p)),2], col="blue",pch=19,bg='grey') # True discoveries
points(total.paths[((total.paths[,3]<=.05) & (total.paths[,2]<bon.p)),5],total.paths[((total.paths[,3]<=.05) & (total.paths[,2]<bon.p)),2], col="green3",pch=19) # Type I errors
points(total.paths[((total.paths[,3]>.05) & (total.paths[,2]>=bon.p)),5],total.paths[((total.paths[,3]>.05) & (total.paths[,2]>=bon.p)),2], col="orange",pch=19) # Type II errors
#

Validated.paths <- total.paths[((total.paths[,3]<=.05) & (total.paths[,2]>=bon.p)),]
save(Validated.paths,file=paste('',dataset,'.GSEA.cv.Rdata',sep="")) # save
save.image(file="Cross_Validation.Rdata")
