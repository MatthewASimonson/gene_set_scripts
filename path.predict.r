gc <- read.table("GIANT.p.05.clumped",header=TRUE)
gc.snps <- as.data.frame(gc$SNP)
names(gc.snps) <- 'SNP'
ga <- read.table("GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt",header=TRUE)
names(ga) <- c('SNP','A1','A2','FRQ','P_g','N')
# simulate genotypes to calculate Beta and standard error for each locus:
# use standard error of 0.02


gc.tot <- merge(gc.snps,ga,by="SNP")
gc.tot <- gc.tot[order(gc.tot$P_g),]
gc.tot$gSTAT <- rep(0,nrow(gc.tot))
gc.tot$gSTAT<- qnorm(gc.tot$P_g/2,lower.tail=FALSE)
gc.tot$BETA <- gc.tot$gSTAT*0.02
profile <- cbind.data.frame(gc.tot$SNP,gc.tot$A1,gc.tot$BETA)
write.table(profile,file="giant.score.05",quote=FALSE,row.names=FALSE,col.names=FALSE)
system("sed -i 's/a/A/g' giant.score.05")
system("sed -i 's/t/T/g' giant.score.05")
system("sed -i 's/c/C/g' giant.score.05")
system("sed -i 's/g/G/g' giant.score.05")
system("plink --bfile REP.SET --score giant.score.05 --pheno BMI.res.phe --out giant.score.05")

#
model1 <- read.table("giant.score.1.profile",header=TRUE)
summary(lm(PHENO~SCORE,data=model1))
#

model05 <- read.table("giant.score.05.profile",header=TRUE)
summary(lm(PHENO~SCORE,data=model05))

model01 <- read.table("giant.score.01.profile",header=TRUE)
summary(lm(PHENO~SCORE,data=model01))

model005 <- read.table("giant.score.005.profile",header=TRUE)
summary(lm(PHENO~SCORE,data=model005))

total <- cbind.data.frame(model1,model05$SCORE,model01$SCORE,model005$SCORE)
names(total) <- c('FID','IID','PHENO','CNT','CNT2','SCORE1','SCORE05','SCORE01','SCORE005') 

summary(lm(PHENO~SCORE1+SCORE05+SCORE01+SCORE005,data=total))
