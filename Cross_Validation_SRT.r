
dataset <- 'MERGE.clean.FINAL'
CORES <- 16
load("resid.models.Rdata") # load total phenotype with covariates for all subjects
total <- total[(!is.na(total$IID) & !is.na(total$BMI) & !is.na(total$SEX) & !is.na(total$AGE) & !is.na(total$BATCH) & !is.na(total$C1)),] # remove any subjects with missing data

  #################################
  # Call multi-threaded packages: #
  #################################

  require(foreach) # 
  require(doMPI) # 
  cl <- startMPIcluster(CORES) 
  registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

  ######################################################################################
  # perform association on random 50% splits of the data in using multi-threaded loop: #
  ######################################################################################
fam <- read.table("MERGE.clean.FINAL.GENE.fam",header=FALSE)
 
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
  
# write phenotype information out for each permutation:
write.table(training.set,file=paste('',dataset,'.train_pheno',i,'',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
write.table(test.set,file=paste('',dataset,'.test_pheno',i,'',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
  #
  
# write out training and test sets for each cross validation:
write.table(training.set[,1:2],file=paste('',dataset,'.train.',i,'',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
write.table(test.set[,1:2],file=paste('',dataset,'.test.',i,'',sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
}
#
  # Permutations:
  loop_train <- foreach(i=1:perms) %dopar% {
    system(paste('plink --bfile ',dataset,'.GENE --keep ',dataset,'.train.',i,' --pheno ',dataset,'.train_pheno',i,' --linear --out ',dataset,'.train',i,' ',sep=""))
    system(paste('plink --bfile ',dataset,'.GENE --keep ',dataset,'.test.',i,' --pheno ',dataset,'.test_pheno',i,' --linear --out ',dataset,'.test',i,' ',sep=""))
  } # END LOOP i
#

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

observed.snp.p <- read.table(paste(dataset,'.assoc.forSRT',sep=""),header=TRUE) 
crossval.data.test <- as.data.frame(matrix(NA,nrow=nrow(observed.snp.p),ncol=perms+2))
crossval.data.train <- as.data.frame(matrix(NA,nrow=nrow(observed.snp.p),ncol=perms+2))

# Store observed associations and SNP name in first 2 columns:
crossval.data.test[,1] <- observed.snp.p$SNP
crossval.data.test[,2] <- observed.snp.p$P
names(crossval.data.test) <- c('SNP','ObservedP')

crossval.data.train[,1] <- observed.snp.p$SNP
crossval.data.train[,2] <- observed.snp.p$P
names(crossval.data.train) <- c('SNP','ObservedP')

# Read in permuted data:
for(h in 1:perms){
print(paste("reading permutation",h,"",sep=""))
crossval.data.train[,h+2] <- read.table(paste(dataset,'.train',h,'.assoc.forSRT',sep=""),header=TRUE,colClasses=c('character','numeric'),nrow=nrow(observed.snp.p))$P
crossval.data.test[,h+2] <- read.table(paste(dataset,'.test',h,'.assoc.forSRT',sep=""),header=TRUE,colClasses=c('character','numeric'),nrow=nrow(observed.snp.p))$P
}

#############################
# Read in pathway with SNP: #
#############################
 
snp.path <- read.table(paste(dataset,'.SRT_paths.txt',sep=""),header=FALSE)
names(snp.path) <- c('PATH','SNP')
path.perms <- merge(snp.path,crossval.data.test,by='SNP') # merge each SNPs p-vals with their pathwway label
path.perms.train <- merge(snp.path,crossval.data.train,by='SNP')
paths <- unique(path.perms$PATH)

########################################################################################################
# Find SNP ratio for each pathway and empirical significance for that pathway BEFORE multiple testing: #
########################################################################################################
 
path.data.test <- as.data.frame(matrix(0,nrow=length(paths),ncol=perms+2))
path.data.train <- as.data.frame(matrix(0,nrow=length(paths),ncol=perms+2))
names(path.data.test) <- c('PATH','Observed.RATIO')
names(path.data.train) <- c('PATH','Observed.RATIO')
threshold <- .01

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
     path.data.test[j,2:ncol(path.data.test)] <- apply(temp.path,2,SR_row,cutoff)

     temp.path.train <- path.perms.train[which(path.perms.train$PATH==paths[j]),3:ncol(path.perms.train)]
     path.data.train[j,2:ncol(path.data.train)] <- apply(temp.path.train,2,SR_row,cutoff)
   } # end loop j
 #
   
path.data.test$PATH <- paths
path.data.train$PATH <- paths

save(path.data.test,file=paste('',dataset,'.path.data.test.Rdata',sep="")) # save object to current directory 
save(path.data.train,file=paste('',dataset,'.path.data.train.Rdata',sep="")) # save object to current directory

# load null from full sample:

load(paste('',dataset,'.path.data.Rdata',sep=""))

# Calculate p-values for each SNP ratio across each cross validation:
 
 get_pval <- function(ob,null){
  pct.a<- mean(null >= ob)
  M <- length(null)
  B1 <- pct.a*M
  pval <- (B1+1)/(M+1)
  return(pval)
  }

train.p.data <- matrix(0,nrow=nrow(path.data.train),ncol=perms+1)
test.p.data <- matrix(0,nrow=nrow(path.data.test),ncol=perms+1)

null.path.data <- cbind.data.frame(path.data[,1],path.data[,3:ncol(path.data)])
#

for(j in 1:perms){
for(i in 1:nrow(null.path.data)){
  SRs <- jitter(as.vector(null.path.data[i,2:ncol(null.path.data)],'numeric'),2)
  SR.ob.train <- path.data.train[i,2+j]
  SR.ob.test <- path.data.test[i,2+j]
  train.p.data[i,j+1] <- get_pval(SR.ob.train,SRs)
  test.p.data[i,j+1] <- get_pval(SR.ob.test,SRs)
  print(paste("pathway",i,""))
} # end loop i
print(paste("validation",j,""))
} # end loop j
#

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
dat <- below_per(train.p.data[,2:ncol(train.p.data)],test.p.data[,2:ncol(test.p.data)],cutoff,path.data[,1])

pval.data <- read.table("MERGE.clean.FINAL.SRT",header=TRUE)
sig.paths <- pval.data[which(pval.data[,2] <= cutoff ),1]
path.names <- path.data[,1]
sig.path.index <- match(sig.paths,path.names)

null.paths <- pval.data[which(pval.data[,2] > cutoff ),1]
null.path.index <- match(null.paths,path.names)
top.paths <- dat[sig.path.index,]
bottom.paths <- dat[null.path.index,]
t.test(top.paths[,2],bottom.paths[,2]) # compare match rates in those above and below cutoff
#

# Create Plot of Data
total.paths.pre <- merge(rbind(top.paths,bottom.paths),pval.data,by='PATHWAY')
tot.pval.index <- order(total.paths.pre[,3])
total.paths <- cbind.data.frame(total.paths.pre[tot.pval.index,],1:512)


plot(total.paths[,2],xlab="Pathway Ranked By Total Sample P-Value",ylab="% Training and Test Both Below Cutoff",main="SRT Pathway Cross-Validation",pch=19,bg='grey')
p <- .05/512
null.05 <- mean(train.p.data[,2:ncol(train.p.data)]<=.05)*mean(test.p.data[,2:ncol(test.p.data)]<=.05) # probability any path cross validation replicates by chance
bon.p <- qexp(p,1/null.05,lower.tail=FALSE) # .05 cutoff before bonferroni


legend(300,.45, # places a legend at the appropriate place
c("Positive","Type I","Type II","Negative","p = 0.05"), # puts text in the legend
pch=c(19,19,19,19,NA),       
lty=c(NA,NA,NA,NA,5), # gives the legend appropriate symbols (lines)
lwd=c(NA,NA,NA,NA,2),       
col=c("blue","green3","orange","black","red"))

points(total.paths[((total.paths[,3]<=.05) & (total.paths[,2]<bon.p)),5],total.paths[((total.paths[,3]<=.05) & (total.paths[,2]<bon.p)),2], col="green3",pch=19) # Type I errors
points(total.paths[((total.paths[,3]>.05) & (total.paths[,2]>=bon.p)),5],total.paths[((total.paths[,3]>.05) & (total.paths[,2]>=bon.p)),2], col="orange",pch=19) # Type II errors
points(total.paths[((total.paths[,3]<=.05) & (total.paths[,2]>=bon.p)),5],total.paths[((total.paths[,3]<=.05) & (total.paths[,2]>=bon.p)),2], col="blue",pch=19,bg='grey') # True discoveries
abline(h=bon.p,col="red",lty=5,lwd=2)
#

# Save pathways that are significant in total sample and pass cross-validation significance 
Validated.paths <- total.paths[((total.paths[,3]<=.05) & (total.paths[,2]>=bon.p)),]
save(Validated.paths,file=paste('',dataset,'.INRICH.cv.Rdata',sep="")) # save object to current directory 

save.image(file="Cross_Validation.Rdata")                  
# Find corrected p-values:
