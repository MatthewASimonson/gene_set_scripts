# Generate INRICH format input files:

dataset <- 'MERGE.clean.FINAL'
PERMS <- 100
CORES <- 8

########################################################################
# 1.) Associated Interval Files for 100 random 50% splits of the data  #
########################################################################
# USE p=.01 cutoff threshold, r^2=.1

  require(foreach) # 
  require(doMPI) # 
  cl <- startMPIcluster(CORES) 
  registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE
  
# Permutations:
  loop_train <- foreach(i=1:PERMS) %dopar% { # for training set  
system(paste("plink --bfile ",dataset," --keep ",dataset,".test.",i," --clump ",dataset,".train",i,".assoc.linear --clump-p1 0.01 --clump-p2 1 --clump-r2 0.10 --clump-kb 250 --clump-range plink.genelist --clump-range-border 20 --out ",dataset,".train",i,"",sep=""))
}

 loop_test <- foreach(i=1:PERMS) %dopar% { # for test set  
system(paste("plink --bfile ",dataset," --keep ",dataset,".train.",i," --clump ",dataset,".test",i,".assoc.linear --clump-p1 0.01 --clump-p2 1 --clump-r2 0.10 --clump-kb 250 --clump-range plink.genelist --clump-range-border 20 --out ",dataset,".test",i,"",sep=""))
}
#

for(i in 1:PERMS){# for training set
# create interval files:
fix.pos <- function(posi){
  posi <- as.character(posi) # convert from factor to character
  loc <- strsplit(posi,':')[[1]][2]
  start.stop <- unlist(strsplit(loc,"\\.."))
  return(start.stop)
} # END FUNCTION fix.pos

# create interval files:
clump <- read.table(paste("",dataset,".train",i,".clumped.ranges",sep=""),header=TRUE)

chrom <- clump$CHR # chromosome information
pos <- as.data.frame(clump$POS)

start.stop <- t(apply(pos,1,fix.pos))

interval <- cbind.data.frame(chrom,start.stop)
write.table(interval,file=paste("",dataset,".train",i,".interval",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}

for(i in 1:PERMS){ # for test set
# create interval files:
fix.pos <- function(posi){
  posi <- as.character(posi) # convert from factor to character
  loc <- strsplit(posi,':')[[1]][2]
  start.stop <- unlist(strsplit(loc,"\\.."))
  return(start.stop)
} # END FUNCTION fix.pos

# create interval files:
clump <- read.table(paste("",dataset,".test",i,".clumped.ranges",sep=""),header=TRUE)

chrom <- clump$CHR # chromosome information
pos <- as.data.frame(clump$POS)

start.stop <- t(apply(pos,1,fix.pos))

interval <- cbind.data.frame(chrom,start.stop)
write.table(interval,file=paste("",dataset,".test",i,".interval",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}

###############################
# 6.) Generate INRICH output  #
 loop_inrich <- foreach(i=1:PERMS) %dopar% { # training set
system(paste("./inrich -g ",dataset,".refgene -m ",dataset,".refSNP -t ",dataset,".gset -a ",dataset,".test",i,".interval -w 20 -r 50000 -q 5000 -i 20 -j 200 -d 0.2 -p 1 -z 1 -o ",dataset,".test",i,"",sep=""))
#
}

 loop_inrich <- foreach(i=1:PERMS) %dopar% { # test set
system(paste("./inrich -g ",dataset,".refgene -m ",dataset,".refSNP -t ",dataset,".gset -a ",dataset,".train",i,".interval -w 20 -r 50000 -q 5000 -i 20 -j 200 -d 0.2 -p 1 -z 1 -o ",dataset,".train",i,"",sep=""))
#
}
#

# Read in results from full sample:

full.inrich <- read.table(paste("",dataset,".out.inrich",sep=""),skip=as.numeric(header.count),nrows=512,fill=TRUE)[,3:5]
names(full.inrich) <- c('Empirical_P','Adjusted_P','PATHWAY')
path.o.index <- order(full.inrich[,1])
pval.data.pre <- full.inrich[path.o.index,]
pval.data.pre[,2] <- gsub("\\*","",as.character(pval.data.pre[,2]))

pval.data <- cbind.data.frame(pval.data.pre[,3],as.numeric(pval.data.pre[,2]),as.numeric(pval.data.pre[,1]))
names(pval.data) <- c('PATHWAY','Adjusted_P','Empirical_P')

# Read in results from cross validations:
paths <- full.inrich$PATHWAY
test.stats <- as.data.frame(matrix(0,nrow=512,ncol=(PERMS))) # matrix to be filled
training.stats <- as.data.frame(matrix(0,nrow=512,ncol=(PERMS))) # matrix to be filled

for (i in 1:PERMS){# read in adjusted p-values for each pathway  
header.count <- 52
train.temp <- read.table(paste("",dataset,".train",i,".out.inrich",sep=""),skip=as.numeric(header.count),nrows=512,fill=TRUE)[,3:5]
test.temp <- read.table(paste("",dataset,".test",i,".out.inrich",sep=""),skip=as.numeric(header.count),nrows=512,fill=TRUE)[,3:5]
training.stats[,i] <- train.temp[,2]
test.stats[,i] <- test.temp[,2]
}
#
train.matrix <- gsub("\\*","",as.matrix(training.stats))
test.matrix <- gsub("\\*","",as.matrix(test.stats))

test.data <- as.data.frame(matrix(0,nrow=512,ncol=(PERMS))) # matrix to be filled
training.data <- as.data.frame(matrix(0,nrow=512,ncol=(PERMS))) # matrix to be filled

for(i in 1:PERMS){
training.data[,i] <- sapply(train.matrix[,i],as.numeric)
test.data[,i] <- sapply(test.matrix[,i],as.numeric)
}

#training.data <- training.stats
#test.data <- test.stats

# Examine the correlation across cross validation samples of p-values from pathways
#########################################################################################

below_per <- function(train,test,cutoff,path.dat){ # returns % of times pathway has p-value below cutoff in both training set AND test set, also p-value for percentage based on empirical null
 train.below <- as.matrix(train) 
 train.below[which(train<=cutoff)] <- 1
 train.below[which(train>cutoff)] <- 0 
 
 test.below <- as.matrix(test)
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
dat <- below_per(training.data,test.data,cutoff,paths)

sig.paths <- pval.data[which(pval.data[,2] <= cutoff ),1]
path.names <- full.inrich$PATHWAY
sig.path.index <- match(sig.paths,path.names)

null.paths <- pval.data[which(pval.data[,2] > cutoff ),1]
null.path.index <- match(null.paths,path.names)
top.paths <- dat[sig.path.index,]
bottom.paths <- dat[null.path.index,]
t.test(top.paths[,2],bottom.paths[,2]) # compare match rates in those above and below cutoff
#

# Create Plot of pathway p-values:

total.paths.pre <- merge(rbind(top.paths,bottom.paths),pval.data,by='PATHWAY')
tot.pval.index <- order(total.paths.pre[,3])
total.paths <- cbind.data.frame(total.paths.pre[tot.pval.index,],1:512)
plot(total.paths[,2],xlab="Pathway Ranked By Total Sample P-Value",ylab="% Training and Test Both Below Cutoff",main="INRICH Pathway Cross-Validation",pch=19,bg='grey')
p <- .05/512
null.05 <- mean(training.data<=.05)*mean(test.data<=.05) # probability any path cross validation replicates by chance
bon.p <- qexp(p,1/null.05,lower.tail=FALSE) # .05 cutoff before bonferroni

legend(300,.45, # places a legend at the appropriate place
c("Positive","Type I","Type II","Negative","p = 0.05"), # puts text in the legend
pch=c(19,19,19,19,NA),       
lty=c(NA,NA,NA,NA,5), # gives the legend appropriate symbols (lines)
lwd=c(NA,NA,NA,NA,2),       
col=c("blue","green3","orange","black","red"))

points(total.paths[((total.paths[,3]<.05) & (total.paths[,2]<bon.p)),5],total.paths[((total.paths[,3]<.05) & (total.paths[,2]<bon.p)),2], col="green3",pch=19) # Type I errors
points(total.paths[((total.paths[,3]>=.05) & (total.paths[,2]>=bon.p)),5],total.paths[((total.paths[,3]>=.05) & (total.paths[,2]>=bon.p)),2], col="orange",pch=19) # Type II errors
points(total.paths[((total.paths[,3]<.05) & (total.paths[,2]>=bon.p)),5],total.paths[((total.paths[,3]<.05) & (total.paths[,2]>=bon.p)),2], col="blue",pch=19,bg='grey') # True discoveries
abline(h=bon.p,col="red",lty=5,lwd=2)

# Save pathways that are significant in total sample and pass cross-validation significance 
Validated.paths <- total.paths[((total.paths[,3]<=.05) & (total.paths[,2]>=bon.p)),]
save(Validated.paths,file=paste('',dataset,'.INRICH.cv.Rdata',sep="")) # save object to current directory 

