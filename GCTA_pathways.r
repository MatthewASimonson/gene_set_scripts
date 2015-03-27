# Generate GRM for each pathway then run GCTA to estimate heritabilities:


paths <- read.table("paths",header=FALSE)

for(i in 1:nrow(paths)){
system(paste("cp h2.sh '",paths$V1[i],".sh'",sep=""))
print(i)
}

for(i in 1:nrow(paths)){
 system(paste("sed -i 's/z/",paths$V1[i],"/g' '",paths$V1[i],".sh'",sep=""))
 print(i)
}

for(i in 1:nrow(paths)){
mrgm.dat1 <- paste("",paths$V1[i],".in",sep="")
mrgm.dat2 <- paste("",paths$V1[i],".out",sep="")
mgrm <- rbind(mrgm.dat1,mrgm.dat2)
write.table(mgrm,file=paste("",paths$V1[i],".mgrm",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
print(i)
}

for(i in 1:nrow(paths)){
system(paste("qsub -q janus-small ",paths$V1[i],".sh",sep=""))
print(i)
} 

for(i in 1:nrow(jid)){
system(paste("qdel ",jid$V1[i],"",sep=""))
print(i)
}

h <- read.table("hsqs",header=FALSE)
VG1 <- vector()
VG1.se <- vector()
VG2 <- vector()
VG2.se <- vector()
SNPs <- vector()
path <- vector()
for(i in 209:nrow(h)){
system(paste("wc -l ",h$V1[i],".snps | cut -d ' ' -f 1 > temp.snp",sep=""))
temp.snp <- read.table("temp.snp",header=FALSE)
SNPs[i] <- temp.snp[1,1]
system(paste("grep 'V(G1)/Vp' ",h$V1[i],".hsq > temp.hsq",sep=""))
temp <- read.table("temp.hsq",header=FALSE)
VG1[i] <- temp[1,2]
VG1.se[i] <- temp[1,3]
system(paste("grep 'V(G2)/Vp' ",h$V1[i],".hsq > temp2.hsq",sep=""))
temp2 <- read.table("temp2.hsq",header=FALSE)
VG2[i] <- temp2[1,2]
VG2.se[i] <- temp2[1,3]
print(i)
path[i] <- as.character(as.matrix(h[i,1])) 
}

V <- cbind.data.frame(VG1,VG1.se,VG2,VG2.se,SNPs,path)
zero.index.p1 <- which(V$VG1.se==0.000000)
zero.index.p2 <- which(V$VG2.se==0.000000)
zero.index <- c(zero.index.p1,zero.index.p2)
V[zero.index,1] <- 0
V[zero.index,2] <- 0

V.nz <- V[-zero.index,]
s.v <- V.nz[order(V.nz$VG1,decreasing=TRUE),] # sort by variance explained
s.v$RATE <- s.v$VG1/s.v$SNPs
s.v <- s.v[order(s.v$VG1,decreasing=TRUE),]
s.v$RANK <- 1:nrow(s.v)

write.table(s.v,file="path.h2.csv",sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

sig.paths <- as.data.frame(c('KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY','KEGG_LYSOSOME','KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY','REACTOME_REGULATION_OF_ORNITHINE_DECARBOXYLASE','REACTOME_STABILIZATION_OF_P53','ST_ERK1_ERK2_MAPK_PATHWAY'))
names(sig.paths) <- c('path')
sig.path.dat <- merge(sig.paths,s.v,by="path")
ns.path.dat <- s.v[-match(sig.path.dat$path,s.v$path),]



# compare heritability explained by significant pathways in discovery sets:

total.snps <- read.table("REP.SET.bim",header=FALSE)
disc.005.sig <- as.matrix(read.table("disc.005.sig",header=FALSE))
disc.01.sig <- as.matrix(read.table("disc.01.sig",header=FALSE))
disc.05.sig <- as.matrix(read.table("disc.05.sig",header=FALSE))
disc.1.sig <- as.matrix(read.table("disc.1.sig",header=FALSE))

#
system("rm disc.005.snps")
for(i in 1:nrow(disc.005.sig)){
system(paste("cat ",disc.005.sig[i,1],".snps >> disc.005.snps",sep=""))
print(i)
}

disc.005.snps <- read.table("disc.005.snps",header=FALSE)
disc.005.u <- unique(disc.005.snps[,1])
write.table(disc.005.u,file="disc.005.in.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)
disc.005.out <- total.snps[-match(disc.005.u,total.snps$V2),2]
write.table(disc.005.out,file="disc.005.out.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)

#
system("rm disc.01.snps")
for(i in 1:nrow(disc.01.sig)){
system(paste("cat ",disc.01.sig[i,1],".snps >> disc.01.snps",sep=""))
print(i)
}

disc.01.snps <- read.table("disc.01.snps",header=FALSE)
disc.01.u <- unique(disc.01.snps[,1])
write.table(disc.01.u,file="disc.01.in.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)
disc.01.out <- total.snps[-match(disc.01.u,total.snps$V2),2]
write.table(disc.01.out,file="disc.01.out.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)

#
system("rm disc.05.snps")
for(i in 1:nrow(disc.05.sig)){
system(paste("cat ",disc.05.sig[i,1],".snps >> disc.05.snps",sep=""))
print(i)
}

disc.05.snps <- read.table("disc.05.snps",header=FALSE)
disc.05.u <- unique(disc.05.snps[,1])
write.table(disc.05.u,file="disc.05.in.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)
disc.05.out <- total.snps[-match(disc.05.u,total.snps$V2),2]
write.table(disc.05.out,file="disc.05.out.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)

#
system("rm disc.1.snps")
for(i in 1:nrow(disc.1.sig)){
system(paste("cat ",disc.1.sig[i,1],".snps >> disc.1.snps",sep=""))
print(i)
}

disc.1.snps <- read.table("disc.1.snps",header=FALSE)
disc.1.u <- unique(disc.1.snps[,1])
write.table(disc.1.u,file="disc.1.in.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)
disc.1.out <- total.snps[-match(disc.1.u,total.snps$V2),2]
write.table(disc.1.out,file="disc.1.out.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Examine null distributin based on 880 pathways:

total.p.sig <- as.matrix(read.table("total.paths",header=FALSE))
total.p.sig <- total.p.sig[-c(218,219),1] # remove non path files
system("rm total.p.snps")
for(i in 1:length(total.p.sig)){
system(paste("cat ",total.p.sig[i]," >> total.p.snps",sep=""))
print(i)
}

total.p.snps <- read.table("total.p.snps",header=FALSE)
total.p.u <- as.character(unique(total.p.snps[,1]))
write.table(total.p.u,file="total.path.in.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)
total.p.out <- total.snps$V2[-match(total.p.u,total.snps$V2)]
write.table(total.p.out,file="path.out.snps",quote=FALSE,row.names=FALSE,col.names=FALSE)

# get size of pathways and genes:

setwd("/STATGEN/home/simonsom/BMI/INRICH/Giant_analysis")
 eg <- read.table("entrez.genes",header=FALSE) # 17838 total genes
 names(eg) <- c('CHR','START','STOP','ENTREZ')
 eg<- eg[-which((eg$CHR=='X')|(eg$CHR=='Y')),] # 16905 autosomal genes

 m <- read.table("msigbd3.0.set",header=TRUE)
 names(m) <- c('ENTREZ','SOURCE','PATH')

 u.genes <- as.data.frame(unique(m$ENTREZ)) # 6804 genes
 names(u.genes) <- c('ENTREZ')
 u.paths <- as.data.frame(unique(m$PATH)) # 880 pathways
 names(u.paths) <- c('PATH')

path.gene.locs <- unique(merge(u.genes,eg,by="ENTREZ"))
tot.length <- sum(path.gene.locs$STOP-path.gene.locs$START)
tot.mb <- tot.length/1000000   # total megabases in pathway genes
avg.path.size <- tot.length/880 # average length of pathway

tp <- rnorm(1000000,mean=7.607388e-05,sd=4.485331e-05) # total pathway variance distribution
# TOTAL PATHWAYS:
# Source  Variance        SE
# V(G1)   0.652801        0.385307
# V(G2)   3.097827        0.793035
# V(e)    19.057555       0.848366
# Vp      22.808183       0.348394
# V(G1)/Vp        0.028621        0.016875
# V(G2)/Vp        0.135821        0.034576
# logL    -17780.954
# n       8621
# # h2 per avg path size: mean = 6.694501e-08 se = 3.947092e-08 


disc.005.paths <- read.table("disc.005.paths",header=FALSE)
names(disc.005.paths) <- c('PATH')
disc.01.paths <- read.table("disc.01.paths",header=FALSE)
names(disc.01.paths) <- c('PATH')
disc.05.paths <- read.table("disc.05.paths",header=FALSE)
names(disc.05.paths) <- c('PATH')
disc.1.paths <- read.table("disc.1.paths",header=FALSE)
names(disc.1.paths) <- c('PATH')

disc.005 <- merge(disc.005.paths,m,by='PATH')
u.005 <- as.data.frame(unique(disc.005$ENTREZ))
names(u.005) <- c('ENTREZ')
p.005.loc <- unique(merge(u.005,eg,by="ENTREZ"))
p.005.length <- sum(p.005.loc$STOP-p.005.loc$START)
scaling.005 <- p.005.length/avg.path.size

p.005 <- rnorm(1000000,mean=9.000912e-05,sd=0.4635876) # total pathway variance distribution
pval.p005 <- 2*pnorm(-1*abs((9.000912e-05)-(6.694501e-08))/(0.4635876)) 
# .005 Pathways:
# Source  Variance        SE
# V(G1)   0.000023        0.118455
# V(G2)   3.787023        0.823975
# V(e)    19.015720       0.848558
# Vp      22.802766       0.348329
# V(G1)/Vp        0.000001        0.005195
# V(G2)/Vp        0.166077        0.035858
# logL    -17780.589
# n       8621
# # h2 per avg path size: mean = 9.000912e-05 se = 0.4635876

disc.01 <- merge(disc.01.paths,m,by='PATH')
u.01 <- as.data.frame(unique(disc.01$ENTREZ))
names(u.01) <- c('ENTREZ')
p.01.loc <- unique(merge(u.01,eg,by="ENTREZ"))
p.01.length <- sum(p.01.loc$STOP-p.01.loc$START)
scaling.01 <- p.01.length/avg.path.size

p.01 <- rnorm(1000000,mean=3.345724e-05,sd=8.769267e-05) # total pathway variance distri
pval.p01 <- 2*pnorm(-1*abs((0.5387438)-(6.694501e-08))/(1.412068)) 
#.01 Pathways:
# Source  Variance        SE
# V(G1)   0.063314        0.165951
# V(G2)   3.682874        0.820238
# V(e)    19.061831       0.848589
# Vp      22.808020       0.348386
# V(G1)/Vp        0.002776        0.007276
# V(G2)/Vp        0.161473        0.035697
# logL    -17780.984
# n       8621
# # h2 per avg path size: mean = 0.5387438 se = 1.412068

disc.05 <- merge(disc.05.paths,m,by='PATH')
u.05 <- as.data.frame(unique(disc.05$ENTREZ))
names(u.05) <- c('ENTREZ')
p.05.loc <- unique(merge(u.05,eg,by="ENTREZ"))
p.05.length <- sum(p.05.loc$STOP-p.05.loc$START)
scaling.05 <- p.05.length/avg.path.size

p.05 <- rnorm(1000000,mean=0.0001379351,sd=8.979071e-05) #
pval.p05 <- 2*pnorm(-1*abs((0.0001379351)-(6.694501e-08))/(8.979071e-05)) 
# .05 Pathways:
# Source  Variance        SE
# V(G1)   0.298229        0.194287
# V(G2)   3.422500        0.820388
# V(e)    19.086430       0.848499
# Vp      22.807158       0.348391
# V(G1)/Vp        0.013076        0.008512
# V(G2)/Vp        0.150063        0.035741
# logL    -17780.559
# n       8621
# # h2 per avg path size: mean = 2.899403 se = 1.887406

disc.1 <- merge(disc.1.paths,m,by='PATH')
u.1 <- as.data.frame(unique(disc.1$ENTREZ))
names(u.1) <- c('ENTREZ')
p.1.loc <- unique(merge(u.1,eg,by="ENTREZ"))
p.1.length <- sum(p.1.loc$STOP-p.1.loc$START)
scaling.1 <- p.1.length/avg.path.size

p.1 <- rnorm(1000000,mean=2.223598,sd=2.012779) #
pval.p1 <- 2*pnorm(-1*abs((2.223598)-(6.694501e-08))/(2.012779)) 
# .1 Pathways:
# Source  Variance        SE
# V(G1)   0.217961        0.197385
# V(G2)   3.541585        0.815886
# V(e)    19.048844       0.848529
# Vp      22.808389       0.348408
# V(G1)/Vp        0.009556        0.008650
# V(G2)/Vp        0.155276        0.035525
# logL    -17780.896
# n       8621
# # h2 per avg path size: mean = 2.223598  se = 2.012779


# Plotting function:


a <-  rnorm(1000000,mean=7.607388e-05,sd=sd(p.005))#/(2*sd(p.005))
e <- p.005#/(2*sd(p.005))
hist(a, 250, prob=TRUE, xlim=c(min(e),max(e)), border="red",col="red",freq=FALSE,axes=TRUE,main=paste("Threshold = Top 0.5%,    p < ",round(pval.p05,3),"",sep=""),xlab="Heritability Per Megabase")
lines(density(a),lwd=5,col='red')
hist(e, 250,prob=TRUE, add=T, border=rgb(0, 0, 1, 0.25),col=rgb(0, 0, 1, 0.5),freq=FALSE,axes=TRUE)
lines(density(e),lwd=3,col='blue')
legend("topright", # places a legend at the appropriate place
c("Expected = 7.6e-05","Observed = 2.64e-08, s.e. = 1.36e-4"), 
fill=c("red","blue"),
bty='n')
#
b <- p.005
c <- p.01
d <- p.05
/sd(c(tp,p.005))


# compare giant genes as a set observed pathways:

 eg <- read.table("entrez.genes",header=FALSE) # 17838 total genes
 names(eg) <- c('CHR','START','STOP','ENTREZ','GENE')
 eg<- eg[-which((eg$CHR=='X')|(eg$CHR=='Y')),1:4] # 16905 autosomal genes

 m <- read.table("msigbd3.0.set",header=TRUE)
 names(m) <- c('ENTREZ','SOURCE','PATH')

 u.genes <- as.data.frame(unique(m$ENTREZ)) # 6804 genes
 names(u.genes) <- c('ENTREZ')
 u.paths <- as.data.frame(unique(m$PATH)) # 880 pathways
 names(u.paths) <- c('PATH')

giant.genes <- read.table("giant.sig.genes",header=FALSE)
names(giant.genes) <- c('ENTREZ','GENE')
gg.2 <- merge(giant.genes,eg,by="ENTREZ")
plink.set <- gg.2[,c(3,4,5,1)]
write.table(plink.set,file="giant.set.genes",quote=FALSE,row.names=FALSE,col.names=FALSE)
system("plink --bfile REP.SET --make-set giant.set.genes --make-set-border 20 --write-set --out giant.genes")
# 1965 SNPs located within giant genes detected with 81 of 90 genes
# Genes without SNPs: 10849,1760,2067,2696,27113,3159,54958,6633,79191

g.in <- read.table("giant.snps.in",header=FALSE)
giant.snps <- read.table("REP.SET.bim",header=FALSE)
giant.snps.out <- giant.snps[-match(g.in$V1,giant.snps$V2),2]
write.table(giant.snps.out,file="giant.snps.out",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("nohup gcta --bfile REP.SET --extract giant.snps.in --make-grm --out giant.gene.in --thread-num 10 &")
system("nohup gcta --bfile REP.SET --extract giant.snps.out --make-grm --out giant.gene.out --thread-num 10 &")

system("nohup gcta --reml --mgrm giant.gene.mgrm --pheno BMI.res.phe --reml-no-lrt --reml-maxit 20 --out giant.genes --thread-num 20 &")

# remove genes without any SNPs
gg.size <- sum(gg.2$STOP-gg.2$START)
