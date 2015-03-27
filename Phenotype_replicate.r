# Phenotype for Replication Set BMI pathway:

# combine information on Gender, Age, BMI, 10 PCAs

# CARDIA and MESA: 
##################

setwd("/home/simonsom/ROH_pathway/TOTAL")
load("phe.model.Rdata")
fam.dat <- read.table("MERGE.clean.FINAL.fam",header=FALSE)
FID.IID <- fam.dat[,1:2]
names(FID.IID) <- c('FID','IID')
# remove NA's from total data
total <- total[(!is.na(total$IID) & !is.na(total$BMI) & !is.na(total$SEX) & !is.na(total$AGE) & !is.na(total$BATCH) & !is.na(total$C1)),]
# remove any duplicates:
dup.index <- duplicated(total$IID)
full.data <- total[which(dup.index==FALSE),]

card.mesa.phe <- full.data[c(which(full.data$set==3),which(full.data$set==2)),]
cmp <- merge(card.mesa.phe,FID.IID,by="IID")

cmp.tot <- cbind.data.frame(cmp$FID,cmp$IID,cmp$SEX,cmp$BATCH,cmp$C1,cmp$C2,cmp$C3,cmp$C4,cmp$C5,cmp$C6,cmp$C7,cmp$C8,cmp$C9,cmp$C10,cmp$BMI,cmp$AGE)
names(cmp.tot) <- c('FID','IID','SEX','PLATE','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','BMI','AGE')

# GNHS: 
#######

setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/GEI")
plate.gei.1 <- read.csv("Sample_annotation_consent_0.csv",header=TRUE)
plate.gei.2 <- read.csv("Sample_annotation_consent_1.csv",header=TRUE)
pca.gei <- read.csv("Principal_components_eur_subset.csv",header=TRUE)

setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/NHS")
plate.nhs <- read.csv("Sample_annotation.csv",header=TRUE)
pca.nhs <- read.csv("Principal_components_white_subset.csv",header=TRUE)

setwd("/home/simonsom/ROH_pathway/TOTAL/GENEVA_diab/TOTAL")
id.key <- read.csv("ID.key.csv",header=TRUE)
ht.g <-read.table("ht.phe",header=TRUE)
bmi.g <-read.csv("GNHS.bmi.total.csv",header=TRUE)
names(ht.g) <- c('FID','SAMPID','HT','AGE')
t1 <- merge(bmi.g,id.key,by="dbGaP.SubjID")
t2 <- merge(t1,ht.g,ny="SAMPID")

ht.g <- cbind.data.frame(t2$FID,t2$SAMPID,t2$bmi,t2$AGE)
names(ht.g) <- c('FID','sample.num','BMI','AGE')
# merge plate and pca from NHS and gei together:

plate.gei <- rbind(plate.gei.1,plate.gei.2)
plate.pca <- merge(plate.gei,pca.gei,by='sample.num')
ppca <- plate.pca[,c(1,2,6,8,22:41)]

plate.pca2 <- merge(plate.nhs,pca.nhs,by='sample.num')
ppca2 <- plate.pca2[,c(1,2,6,8,26:45)]
names(ppca2) <- names(ppca)

ppca.tot <- rbind(ppca,ppca2)

geneva.dat <- merge(ppca.tot,ht.g,by='sample.num')
geneva.dat <- cbind.data.frame(geneva.dat,rep(6,nrow(geneva.dat)))
names(geneva.dat) <- c('IID','GID','SEX','PLATE','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','FID','BMI','AGE','set')

GNHS.dat <- geneva.dat[,c(1,3:14,26,27)]
GNHS.dat <- cbind.data.frame(rep(-1,nrow(GNHS.dat)),GNHS.dat)
names(GNHS.dat) <- c('FID','IID','SEX','PLATE','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','BMI','AGE')
GNHS.dat <- tf[,c(2,1,8:ncol(tf))]
# MEGE all phenotype data from all sets:

gmc <- rbind(as.matrix(cmp.tot),as.matrix(GNHS.dat)) 
gmc <- as.data.frame(gmc[complete.cases(gmc),])
setwd("/home/simonsom/BMI/INRICH/Replicate")

write.table(gmc,file="Phenotype",quote=FALSE,row.names=FALSE,col.names=TRUE)

# Residualize out covariates:

bmi.model <- lm(as.numeric(BMI) ~ as.factor(SEX) + as.factor(PLATE) + as.numeric(AGE) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) + as.numeric(C6) + as.numeric(C7) + as.numeric(C8) + as.numeric(C9) + as.numeric(C10),data=gmc)

bmi.res.phe <- cbind.data.frame(gmc$FID,gmc$IID,resid(bmi.model))
names(bmi.res.phe) <- c('FID','IID','BMI')
write.table(bmi.res.phe,file="BMI.res.phe",quote=FALSE,row.names=FALSE,col.names=TRUE)
