# Get interval index SNPs from each pathway's intervals then run fixed effects model on each pathway:



                    
# Calculate h2 SNP estimates for total BMI:

tot.path.snps.pre <- read.table("path.snps",header=FALSE)
tot.path.snps <- unique(tot.path.snps.pre)
total.snps <- read.table("REP.SET.bim",header=FALSE)
path.snp.index <- match(total.snps$V2,tot.path.snps$V1)
np.snp.index <- which(is.na(path.snp.index)==TRUE)
non.path.snps <- total.snps[np.snp.index,2]
write.table(non.path.snps,file="non.path.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)
#

system(paste("nohup gcta --bfile REP.SET --extract path.snps --make-grm --out path.h2.snp &",sep=""))
system(paste("nohup gcta --bfile REP.SET --extract non.path.snps --make-grm --out non.path.h2.snp &",sep=""))






