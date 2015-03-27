# Generate INRICH format input files:

dataset <- 'MERGE.clean.FINAL'

############################
# 1.) Target Gene Set File #
############################

 load(paste('',dataset,'.gene.set.Rdata',sep="")) # read in gene set information created using 'generate_gene_set()' function

 gene.set <- gene.sets.final

gene_set_frame <- function(gene.set){
 frame.rows <- sum(unlist(lapply(gene.set,nrow))) # number of unique gene/pathway combinations

 # fill path column of matrix:
  gene.path.mat <- matrix('X',nrow=frame.rows,ncol=2)
  path.count <- lapply(gene.set,nrow)
 
  for(i in 1:length(path.count)){
     start <- sum(as.numeric(path.count[1:i-1]))+1
     end <- (start-1)+as.numeric(path.count[[i]])
     gene.path.mat[start:end,2] <- names(path.count[i])
     gene.path.mat[start:end,1] <- as.character(unlist(gene.set[[i]]))
  }# END LOOP i

frame <- as.data.frame(gene.path.mat) 
 return(frame)
} # END FUNCTION gene_set_frame

frame <- gene_set_frame(gene.set)
frame$V3 <- rep('PATH ANNOTATION NA',nrow(frame))

write.table(frame,file=paste('',dataset,'.gset',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

###########################
# 2.) Reference Gene File #
###########################

pre.ref <- read.table("clean.genelist.FINAL",header=TRUE) # read in output from function 'clean_gene_list()'

ref.gene <- cbind.data.frame(pre.ref$CHR,pre.ref$START,pre.ref$STOP,pre.ref$NAME,pre.ref$NAME)
ref.gene$ANNOT <- rep('PATH_ANNOTATION_NA',nrow(ref.gene))

write.table(ref.gene,file=paste('',dataset,'.refgene',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

##########################
# 3.) Reference SNP File #
##########################

pre.SNP.ref <- read.table(paste('',dataset,'.bim',sep=""),header=FALSE)
SNP.ref <- cbind.data.frame(pre.SNP.ref$V1,pre.SNP.ref$V4)
write.table(SNP.ref,file=paste('',dataset,'.refSNP',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
SNPlist <- as.data.frame(pre.SNP.ref$V2)
write.table(SNPlist,file=paste('',dataset,'.SNPlist',sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

################################
# 4.) Gene List file           #
################################

pre.gen <- read.table("clean.genelist.FINAL",header=TRUE) # read in output from function 'clean_gene_list()'
plink.genelist <- cbind.data.frame(pre.gen[,2:4],pre.gen[,1])
write.table(plink.genelist,file='plink.genelist',row.names=FALSE,col.names=FALSE,quote=FALSE)

###############################################
# 5.) Associated Interval Files for real data #
###############################################
# USE p=.01 cutoff threshold, r^2=.1
system(paste("plink --bfile ",dataset," --clump ",dataset,".assoc.linear --clump-p1 0.01 --clump-p2 1 --clump-r2 0.10 --clump-kb 250 --clump-range plink.genelist --clump-range-border 20 --out ",dataset,"",sep=""))

# create interval files:
fix.pos <- function(posi){
  posi <- as.character(posi) # convert from factor to character
  loc <- strsplit(posi,':')[[1]][2]
  start.stop <- unlist(strsplit(loc,"\\.."))
  return(start.stop)
} # END FUNCTION fix.pos

# create interval files:
clump <- read.table(paste("",dataset,".clumped.ranges",sep=""),header=TRUE)

chrom <- clump$CHR # chromosome information
pos <- as.data.frame(clump$POS)

start.stop <- t(apply(pos,1,fix.pos))

interval <- cbind.data.frame(chrom,start.stop)
write.table(interval,file=paste("",dataset,".interval",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)

###############################
# 6.) Generate INRICH output  #

system(paste("./inrich -g ",dataset,".refgene -m ",dataset,".refSNP -t ",dataset,".gset -a ",dataset,".interval -w 20 -r 50000 -q 5000 -i 20 -j 200 -d 0.2 -p 1 -z 1 -o ",dataset,"",sep=""))
#

