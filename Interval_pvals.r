# custom regional association plot for intervals


custom_ra_plot<- function(locus, map, genes, int.start, int.stop, flanking = 1000, best.pval = NULL, sf = c(4, 4), logpmax = 10, pch = 21){
    NAME <- locus$NAME
    PVAL <- locus$PVAL
    snp <- NAME[PVAL == min(PVAL[which((locus$POS>=int.start)&(locus$POS<=int.stop))])] # min pval inside interval
    chr <- locus$CHR[1]
    hit <- locus[NAME == snp, ]
    if (is.null(best.pval)){ 
        best.pval <- hit$PVAL
      }
    lb <- min(locus$POS) - flanking
    ub <- max(locus$POS) + flanking
    lu <- ub - lb
    center <- lb + lu/2
    center.100kb <- round(center/1e+05) * 1e+05
    offset.100kb <- round(lu/5/1e+05) * 1e+05
    offset <- logpmax/sf[1]
    ylim <- logpmax + offset
    ylim4 <- ylim/4
    yadj <- -offset + ylim/sf[2]
    keep <- subset(map, map$POS > lb & map$POS < ub)
    genes <- with(genes, subset(genes, (START > lb & START < 
        ub) | (STOP > lb & STOP < ub)))
    par(mar = c(4, 4, 3, 4))
    xy <- xy.coords(keep$POS, yadj + (keep$THETA/60) * (3 * ylim4))
    plot(xy$x, xy$y, type = "l", col = "darkblue", lwd = 1, 
        ylim = c(-offset, logpmax), xlab = "", ylab = "", axes = F)
    box()
    if (offset.100kb != 0) 
        center5 <- center.100kb + (-2:2) * offset.100kb
    else {
        p1 <- min(xy$x)
        p5 <- max(xy$x)
        p3 <- (p1 + p5)/2
        p2 <- (p1 + p3)/2
        p4 <- (p3 + p5)/2
        center5 <- c(p1, p2, p3, p4, p5)
    }
    axis(1, at = center5, labels = round(center5/1000, 2), las = 1)
    mtext(paste("Chromosome", chr, "position (kb)", sep = " "), 
        side = 1, line = 2.5)
    axis(2, at = seq(0, logpmax, 2), labels = seq(0, logpmax, 
        2), las = 1)
    mtext("-log10(Observed p)", side = 2, at = logpmax/2, line = 2)
    axis(4, at = yadj + (0:4) * ylim4, labels = paste((0:4) * 
        20), las = 1, col.axis = "darkblue")
    mtext("Recombination rate (cM/Mb)", side = 4, at = -offset + 
        2 * ylim4, line = 2, col = "darkblue")
    lines(c(lb, ub), c(0, 0), lty = "dotted", lwd = 1, col = "black")
    lines(c(int.start, int.stop), c(0, 0), lty = "solid", lwd = 3, col = "deeppink2")
    points(hit$POS, -log10(hit$PVAL), pch = pch, cex = 2.5, bg = "red")
    text(hit$POS, -log10(hit$PVAL), labels = hit$NAME, pos = 3, 
        offset = 1)
    if (-log10(best.pval) < logpmax) {
        points(hit$POS, -log10(best.pval), pch = pch, cex = 2.5, 
            bg = "blue")
        text(hit$POS, -log10(best.pval), labels = c(paste("P=", 
            best.pval, sep = "")), pos = 4, offset = 2)
    }
    else {
        points(hit$POS, logpmax, pch = pch, cex = 2.5, bg = "blueviolet")
        text(hit$POS, logpmax, labels = c(paste("P=", best.pval, 
            sep = "")), pos = 4, offset = 1)
    }
    RSQR <- locus$RSQR
    strong.ld <- subset(locus, NAME != snp & RSQR >= 0.8)
    moderate.ld <- subset(locus, NAME != snp & RSQR >= 0.5 & 
        RSQR < 0.8)
    weak.ld <- subset(locus, NAME != snp & RSQR >= 0.2 & RSQR < 
        0.5)
    not.in.ld <- subset(locus, NAME != snp & RSQR < 0.2)
    na.ld <- subset(locus, NAME != snp & is.na(RSQR))
    colors <- c("red","orange","yellow","grey","white")
    points(strong.ld$POS, -log10(strong.ld$PVAL), pch = pch, 
        cex = 1.25, bg = colors[1])
    points(moderate.ld$POS, -log10(moderate.ld$PVAL), pch = pch, 
        cex = 1.25, bg = colors[2])
    points(weak.ld$POS, -log10(weak.ld$PVAL), pch = pch, cex = 1.25, 
        bg = colors[3])
    points(not.in.ld$POS, -log10(not.in.ld$PVAL), pch = pch, 
        cex = 1, bg = colors[4])
    points(na.ld$POS, -log10(na.ld$PVAL), pch = pch, cex = 1, 
        bg = colors[5])
    for (i in 1:nrow(genes)) {
        gi <- genes[i, ]
        GENE <- gi$GENE
        cat("-", GENE, "\n")
        start <- gi$START
        finish <- gi$STOP
        center <- (start + finish)/2
        lb <- min(xy$x)
        ub <- max(xy$x)
        adj <- -offset + 2 * (i%%4)/3
        if (!is.na(GENE)) {
            if (start < lb) 
                text(lb, adj + ylim4/10, labels = GENE, cex = 0.7)
            else if (finish > ub) 
                text(ub, adj + ylim4/10, labels = GENE, cex = 0.7)
            else text(center, adj + ylim4/10, labels = GENE, 
                cex = 0.7)
        }
        STRAND <- gi$STRAND
        if (STRAND == "+") 
            arrows(start, adj, finish, adj, length = 0.05, lwd = 2, 
                code = 2, lty = "solid", col = "darkgreen")
        else arrows(start, adj, finish, adj, length = 0.05, lwd = 2, 
            code = 1, lty = "solid", col = "darkgreen")
    }
    ltext <- rbind("0.8-", "0.5-", "0.2-", "0.0-", "NA")
    legend(min(xy$x), logpmax, legend = ltext, title = expression(r^2), 
        fill = colors)
}

###############################################
# Examine intervals in each significant pathway

# Interval p-values:
system("grep '_O2' REP.p.005.out.inrich | cut -f 2-4 > REP.p.005.int.stats1")
system("grep '_O2' REP.p.005.out.inrich | cut -f 6- | cut -d ' ' -f 2- > REP.p.005.int.stats2")
system("paste REP.p.005.int.stats1 REP.p.005.int.stats2 > REP.p.005.int.stats")
r <- read.table("REP.p.005.int.stats",header=TRUE)
s <- read.table("REP.p.005.clumped.ranges",header=TRUE)
s$INTERVAL <- paste("int",0:(nrow(s)-1),":",s$POS,"",sep="")
inrich.int.stats <- merge(r,s,by="INTERVAL")
in.stats <- inrich.int.stats[order(inrich.int.stats$TARGET,inrich.int.stats$GENE_LOC),] # summary data for all intervals in rep paths .005

gwas <- read.table("REP.GWAS.qassoc",header=TRUE) # read in GWAS from replication set
#gwas.snps <- gwas[,c(2,9)]
#write.table(gwas.snps,file="REP.gwas.csv",quote=FALSE,row.names=FALSE,col.names=TRUE,sep=",")

kli <- in.stats[in.stats$TARGET=='KEGG_LYSOSOME',]

Feature <- vector()
chr <- vector()
start <- vector()
end <- vector()
flank <- vector()
plot <- vector()
arguments <- vector()

#system("cut -d ' ' -f 4,5 entrez.gene.map > entrez.key") 
#system("cut -d ' ' -f 1,2,3,5 entrez.gene.map > entrez.gene.locs")

entrez.key <- read.table("entrez.key",header=FALSE)
kli.hugo <- as.character(entrez.key[match(kli$GENE_ID,entrez.key[,1]),2]) 
  
for(i in 1:nrow(kli)){
Feature[i] <- kli.hugo[i]
chr[i] <- as.numeric(kli$CHR[i])
end[i] <- as.numeric(strsplit(as.character(kli$POS[i]),"\\..")[[1]][2])#as.numeric(strsplit(as.character(kli$GENE_LOC[i]),"\\..")[[1]][2])
front <- strsplit(as.character(kli$POS[i]),"\\..")[[1]][1] #strsplit(as.character(kli$GENE_LOC[i]),"\\..")[[1]][1]
start[i] <- as.numeric(strsplit(front,":")[[1]][2])#as.numeric(strsplit(front,":")[[1]][2])
i.end <- as.numeric(strsplit(as.character(kli$POS[i]),"\\..")[[1]][2])
i.front <- strsplit(as.character(kli$POS[i]),"\\..")[[1]][1]
i.start <- as.numeric(strsplit(i.front,":")[[1]][2])
flank[i] <- '250kb'
plot[i] <- 'yes'
arguments[i] <- paste('',kli.hugo[i],' Interval',sep='')
}

kli.hitspec <- cbind.data.frame(Feature,chr,start,end,flank,plot,arguments,kli$SNP) 

# Calculate LD information for SNPs surrounding index:
for(i in 1:nrow(kli)){
system(paste("nohup plink --bfile REP.SET --r2 --ld-snp ",kli$SNP[i]," --ld-window-kb 2000 --ld-window 99999 --ld-window-r2 0 --out ",kli$SNP[i]," &",sep=""))
}
#
gene.locs <- read.table("entrez.gene.locs",header=FALSE)
p.thresh <- .005 # set pval thresh here
# create regional association plots for each interval:
for(j in 3:nrow(kli.hitspec)){
ld <- read.table(paste("",kli$SNP[j],".ld",sep=""),header=TRUE)
rel.ld <- ld[,c(6,7)]
names(rel.ld) <- c('SNP','RSQR')


names(gene.locs) <- c('CHR','START','STOP','GENE')
genes <- gene.locs[gene.locs$CHR==kli.hitspec$chr[j],]
gene.dat <- cbind.data.frame(genes$START,genes$STOP,rep("-",nrow(genes)),genes$GENE)
names(gene.dat) <- c('START','STOP','STRAND','GENE')

chr.index <- which(gwas$CHR==kli.hitspec$chr[j])
locus <- gwas[chr.index,c(1,3,2,9)]
locus <- merge(locus,rel.ld,by="SNP")
locus <- locus[,c(2,3,1,4,5)]
names(locus) <- c('CHR','POS','NAME','PVAL','RSQR')
loc.keep <- which(locus$POS>=(kli.hitspec$start[j]-50000) & (locus$POS<=(kli.hitspec$end[j]+50000)))
locus.dat <- locus[loc.keep,] # locus data surrounding relevant SNP
min.in <- min(locus.dat$POS[locus.dat$POS>=(kli.hitspec$start[j])])
if(locus.dat[which(locus.dat$POS==min.in),5]<.2){
 locus.dat[which(locus.dat$POS==min.in),5] <- .2
}
max.in <- max(locus.dat$POS[locus.dat$POS<=(kli.hitspec$end[j])])
if(locus.dat[which(locus.dat$POS==max.in),5]<.2){
 locus.dat[which(locus.dat$POS==max.in),5] <- .2
}
low.out <- which(locus.dat$POS<kli.hitspec$start[j])
high.out <- which(locus.dat$POS>kli.hitspec$end[j])
locus.dat[c(low.out,high.out),5] <- .1

map <- read.table(paste("genetic_map_chr",kli.hitspec$chr[j],"_b36.txt",sep=""),header=TRUE)
names(map) <- c('POS','THETA','DIST')
map.keep <- which((map$POS>=(kli.hitspec$start[j]-50000) & (map$POS<=kli.hitspec$end[j])+50000))
map.dat <- map[map.keep,] 

pdf(file=paste("",kli.hitspec$Feature[j],".pct",as.character(p.thresh),".pdf",sep=""))
custom_ra_plot(locus.dat, map.dat, gene.dat,kli.hitspec$start[j],kli.hitspec$end[j], flanking=1000, best.pval=locus.dat$PVAL[which(as.character(locus.dat$NAME)==as.character(kli$SNP[j]))], sf=c(4,4), logpmax=10, pch=21)
title(paste("Interval ",kli.hitspec$Feature[j]," using cutoff of top ",as.character(p.thresh),"%",sep=""))
#
dev.off()
print(j)
}
#

require(foreach) # 
require(doMPI) # 
cl <- startMPIcluster(CORES) 
registerDoMPI(cl) # DoMC must be registered so the foreach tasks are executed in parallel instead of sequentially; NUMBER OF CORES CAN BE SET HERE

loop_GRM <- foreach(i=1:nrow(bin.genes)) %dopar% {
  gene <- as.character(bin.genes$NAME[i])
  system(paste("plink --bfile ",gene," --noweb --recode12 --out ",gene,"",sep="")) # recode data as non-binary and remove subjects without SNPs
 } # end loop_GENERATE_PED
# j is 8, still has to be completed
