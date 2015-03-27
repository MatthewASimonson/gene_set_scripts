########################################################
# Generate file with SNPs in linkage equilibrian for SRT
########################################################

dataset <- 'ARIC.clean'

system(paste("plink --bfile ",dataset," --indep 50 5 2 --out ",dataset,"",sep="")) 

system(paste("plink --bfile ",dataset," --extract ",dataset,".prune.in --make-bed --out ",dataset,".LE",sep=""))
