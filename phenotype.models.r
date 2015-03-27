# Workspace image saved in working directory with all variables

load("phe.model.Rdata")
fam.dat <- read.table("MERGE.clean.FINAL.fam",header=FALSE)
FID.IID <- fam.dat[,1:2]
names(FID.IID) <- c('FID','IID')
# remove NA's from total data
total <- total[(!is.na(total$IID) & !is.na(total$BMI) & !is.na(total$SEX) & !is.na(total$AGE) & !is.na(total$BATCH) & !is.na(total$C1)),]
# remove any duplicates:
dup.index <- duplicated(total$IID)
full.data <- total[which(dup.index==FALSE),]

# total data model:
total.mod <- lm(as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) +as.numeric(C6)+ as.numeric(C7)+ as.numeric(C8)+ as.numeric(C9)+ as.numeric(C10)+ as.numeric(C11)+ as.numeric(C12)+ as.numeric(C13)+ as.numeric(C14)+ as.numeric(C15)+ as.numeric(C16) +as.numeric(C17)+ as.numeric(C18)+as.numeric(C19) +as.numeric(C20)+as.factor(set)+as.factor(BATCH),data=full.data)

rt <- resid(total.mod)
krt <- kurtosis(rt)
skw <- skewness(rt)
phe.total <- cbind.data.frame(full.data$IID,rt)
names(phe.total) <- c('IID','BMI')
phe.t <- merge(FID.IID,phe.total,by="IID")
write.table(phe.t[,c(2,1,3)],file="MERGE.clean.FINAL.phe",quote=FALSE,row.names=FALSE,col.names=FALSE)

###############################################
# For ARIC:

# remove NA's from total data
a <- a[(!is.na(a$IID) & !is.na(a$BMI) & !is.na(a$SEX) & !is.na(a$AGE) & !is.na(a$BATCH) & !is.na(a$C1)),]

# ARIC data model:
ARIC.mod <- lm(as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) +as.numeric(C6)+ as.numeric(C7)+ as.numeric(C8)+ as.numeric(C9)+ as.numeric(C10)+ as.numeric(C11)+ as.numeric(C12)+ as.numeric(C13)+ as.numeric(C14)+ as.numeric(C15)+ as.numeric(C16) +as.numeric(C17)+ as.numeric(C18)+as.numeric(C19) +as.numeric(C20)+as.factor(BATCH),data=a)

rt.a <- resid(ARIC.mod)
phe.total.a <- cbind.data.frame(a$IID,rt.a)
names(phe.total.a) <- c('IID','BMI')
phe.ta <- merge(FID.IID,phe.total.a,by="IID")
write.table(phe.ta[,c(2,1,3)],file="ARIC.clean.FINAL.phe",quote=FALSE,row.names=FALSE,col.names=FALSE)

###############################################
# For CARDIA:

# remove NA's from total data
c$SEX[which(c$SEX==2)] <- 'F'
c$SEX[which(c$SEX==1)] <- 'M'

c <- c[(!is.na(c$IID) & !is.na(c$BMI) & !is.na(c$SEX) & !is.na(c$AGE) & !is.na(c$BATCH) & !is.na(c$C1)),]

# CARDIA data model:
CARDIA.mod <- lm(as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) +as.numeric(C6)+ as.numeric(C7)+ as.numeric(C8)+ as.numeric(C9)+ as.numeric(C10)+ as.numeric(C11)+ as.numeric(C12)+ as.numeric(C13)+ as.numeric(C14)+ as.numeric(C15)+ as.numeric(C16) +as.numeric(C17)+ as.numeric(C18)+as.numeric(C19) +as.numeric(C20)+as.factor(BATCH),data=c)

rt.c <- resid(CARDIA.mod)
phe.total.c <- cbind.data.frame(c$IID,rt.c)
names(phe.total.c) <- c('IID','BMI')
phe.tc <- merge(FID.IID,phe.total.c,by="IID")
write.table(phe.tc[,c(2,1,3)],file="CARDIA.clean.FINAL.phe",quote=FALSE,row.names=FALSE,col.names=FALSE)
###############################################
# For MESA:
# remove NA's from total data
m$SEX[which(m$SEX==2)] <- 'F'
m$SEX[which(m$SEX==1)] <- 'M'
# remove NA's from total data
m <- m[(!is.na(m$IID) & !is.na(m$BMI) & !is.na(m$SEX) & !is.na(m$AGE) & !is.na(m$BATCH) & !is.na(m$C1)),]

# MESA data model:
MESA.mod <- lm(as.numeric(BMI) ~ as.numeric(AGE) + as.factor(SEX) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) +as.numeric(C6)+ as.numeric(C7)+ as.numeric(C8)+ as.numeric(C9)+ as.numeric(C10)+ as.numeric(C11)+ as.numeric(C12)+ as.numeric(C13)+ as.numeric(C14)+ as.numeric(C15)+ as.numeric(C16) +as.numeric(C17)+ as.numeric(C18)+as.numeric(C19) +as.numeric(C20)+as.factor(BATCH),data=m)

rt.m <- resid(MESA.mod)
phe.total.m <- cbind.data.frame(m$IID,rt.m)
names(phe.total.m) <- c('IID','BMI')
phe.tm <- merge(FID.IID,phe.total.m,by="IID")
write.table(phe.tm[,c(2,1,3)],file="MESA.clean.FINAL.phe",quote=FALSE,row.names=FALSE,col.names=FALSE)
