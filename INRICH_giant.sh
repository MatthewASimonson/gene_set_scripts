# Generate INRICH input files:
################################
# Extract pathways that contain greater than 1 interval below threshold

for p005 in `cat top.005.paths`
do
        grep -w $p005 msigbd3.0.set >> p005.paths
done

for p01 in `cat top.01.paths`
do
        grep -w $p01 msigbd3.0.set >> p01.paths
done

for p05 in `cat top.05.paths`
do
        grep -w $p05 msigbd3.0.set >> p05.paths
done

for p1 in `cat top.1.paths`
do
        grep -w $p1 msigbd3.0.set >> p1.paths
done
# create reference SNP file:

gawk ' NR>1 { print $1,$4 }' hapmap-ceu-commons.bim  > GIANT.map 

# Generate interval file from hapmap: (only includes SNPs that are common across all data sets)
# Cite Lee Depression paper for parameters
plink --bfile hapmap-ceu-commons \
      --clump GIANT.assoc \
      --clump-p1 0.00137 \
      --clump-p2 1 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --clump-range entrez.gene.map \
      --clump-range-border 20 \
      --out GIANT.p.005
gawk 'NR>1 {print $1 "." substr($5, index($5,":")+1) }' \
GIANT.p.005.clumped.ranges | gawk -F. '{ print $1,$2,$4 }' > GIANT.p.005.interval

# Run INRICH using all SNPs on msigdb 3.0 Pathways on GIANT discovery set:

./inrich \
-a GIANT.p.005.interval \
-m GIANT.map \
-g entrez.gene.map \
-t p005.paths \
-w 20 \
-r 100000 \
-q 500 \
-i 20 \
-j 200 \
-p 1 \
-d .1 \
-z 2 \
-o GIANT.p.005

# Generate interval file from hapmap: (only includes SNPs that are common across all data sets)
# Cite Lee Depression paper for parameters
plink --bfile hapmap-ceu-commons \
      --clump GIANT.assoc \
      --clump-p1 0.00432 \
      --clump-p2 1 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --clump-range entrez.gene.map \
      --clump-range-border 20 \
      --out GIANT.p.01
gawk 'NR>1 {print $1 "." substr($5, index($5,":")+1) }' \
GIANT.p.01.clumped.ranges | gawk -F. '{ print $1,$2,$4 }' > GIANT.p.01.interval

# Run INRICH using all SNPs on msigdb 3.0 Pathways on GIANT discovery set:

./inrich \
-a GIANT.p.01.interval \
-m GIANT.map \
-g entrez.gene.map \
-t p01.paths \
-w 20 \
-r 100000 \
-q 5000 \
-i 20 \
-j 200 \
-p 1 \
-d .1 \
-z 2 \
-o GIANT.p.01

# Generate interval file from hapmap: (only includes SNPs that are common across all data sets)
# Cite Lee Depression paper for parameters
plink --bfile hapmap-ceu-commons \
      --clump GIANT.assoc \
      --clump-p1 0.0385 \
      --clump-p2 1 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --clump-range entrez.gene.map \
      --clump-range-border 20 \
      --out GIANT.p.05
gawk 'NR>1 {print $1 "." substr($5, index($5,":")+1) }' \
GIANT.p.05.clumped.ranges | gawk -F. '{ print $1,$2,$4 }' > GIANT.p.05.interval

# Run INRICH using all SNPs on msigdb 3.0 Pathways on GIANT discovery set:

nohup ./inrich \
-a GIANT.p.05.interval \
-m GIANT.map \
-g entrez.gene.map \
-t p05.paths \
-w 20 \
-r 100000 \
-q 5000 \
-i 20 \
-j 200 \
-p 1 \
-d .1 \
-z 2 \
-o GIANT.p.05 &

# Generate interval file from hapmap: (only includes SNPs that are common across all data sets)
# Cite Lee Depression paper for parameters
plink --bfile hapmap-ceu-commons \
      --clump GIANT.assoc \
      --clump-p1 0.087 \
      --clump-p2 1 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --clump-range entrez.gene.map \
      --clump-range-border 20 \
      --out GIANT.p.1
gawk 'NR>1 {print $1 "." substr($5, index($5,":")+1) }' \
GIANT.p.1.clumped.ranges | gawk -F. '{ print $1,$2,$4 }' > GIANT.p.1.interval

# Run INRICH using all SNPs on msigdb 3.0 Pathways on GIANT discovery set:

nohup ./inrich \
-a GIANT.p.1.interval \
-m GIANT.map \
-g entrez.gene.map \
-t msigbd3.0.set \
-w 20 \
-r 100000 \
-q 5000 \
-i 20 \
-j 200 \
-p 1 \
-d .1 \
-z 2 \
-o GIANT.p.1 &

###########################
# Replication Set INRICH: #
###########################

# Generate interval file from hapmap: (only includes SNPs that are common across all data sets)
# Cite Lee Depression paper for parameters
plink --bfile hapmap-ceu-commons \
      --clump REP.GWAS.qassoc \
      --clump-p1 0.004272 \
      --clump-p2 1 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --clump-range entrez.gene.map \
      --clump-range-border 20 \
      --out REP.p.005
gawk 'NR>1 {print $1 "." substr($5, index($5,":")+1) }' \
REP.p.005.clumped.ranges | gawk -F. '{ print $1,$2,$4 }' > REP.p.005.interval

# Run INRICH using all SNPs on msigdb 3.0 Pathways on REP set:
for rep005 in `cat top.005.paths.rep`
do
        grep -w $rep005 msigbd3.0.set >> rep005.paths
done

./inrich \
-a REP.p.005.interval \
-m GIANT.map \
-g entrez.gene.map \
-t rep005.paths \
-w 20 \
-r 100000 \
-q 5000 \
-i 20 \
-j 200 \
-p 1 \
-d .1 \
-o REP.p.005

# Generate interval file from hapmap: (only includes SNPs that are common across all data sets)
# Cite Lee Depression paper for parameters
plink --bfile hapmap-ceu-commons \
      --clump REP.GWAS.qassoc \
      --clump-p1 0.008656 \
      --clump-p2 1 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --clump-range entrez.gene.map \
      --clump-range-border 20 \
      --out REP.p.01
gawk 'NR>1 {print $1 "." substr($5, index($5,":")+1) }' \
REP.p.01.clumped.ranges | gawk -F. '{ print $1,$2,$4 }' > REP.p.01.interval

# Run INRICH using all SNPs on msigdb 3.0 Pathways on REP set:
rm p01.paths
for rep01 in `cat top.01.paths.pre`
do
        grep -w $rep01 msigbd3.0.set >> p01.paths
done

./inrich \
-a REP.p.01.interval \
-m GIANT.map \
-g entrez.gene.map \
-t p01.paths \
-w 20 \
-r 100000 \
-q 5000 \
-i 20 \
-j 200 \
-p 1 \
-d .1 \
-o REP.p.01


# Generate interval file from hapmap: (only includes SNPs that are common across all data sets)
# Cite Lee Depression paper for parameters
plink --bfile hapmap-ceu-commons \
      --clump REP.GWAS.qassoc \
      --clump-p1 0.0456 \
      --clump-p2 1 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --clump-range entrez.gene.map \
      --clump-range-border 20 \
      --out REP.p.05
gawk 'NR>1 {print $1 "." substr($5, index($5,":")+1) }' \
REP.p.05.clumped.ranges | gawk -F. '{ print $1,$2,$4 }' > REP.p.05.interval

# Run INRICH using all SNPs on msigdb 3.0 Pathways on REP set:
rm p05.paths
for rep1 in `cat top.05.paths.pre`
do
        grep -w $rep1 msigbd3.0.set >> p05.paths
done

nohup ./inrich \
-a REP.p.05.interval \
-m GIANT.map \
-g entrez.gene.map \
-t p05.paths \
-w 20 \
-r 100000 \
-q 5000 \
-i 20 \
-j 200 \
-p 1 \
-d .1 \
-o REP.p.05 &

# Generate interval file from hapmap: (only includes SNPs that are common across all data sets)
# Cite Lee Depression paper for parameters
plink --bfile hapmap-ceu-commons \
      --clump REP.GWAS.qassoc \
      --clump-p1 0.0938 \
      --clump-p2 1 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --clump-range entrez.gene.map \
      --clump-range-border 20 \
      --out REP.p.1
gawk 'NR>1 {print $1 "." substr($5, index($5,":")+1) }' \
REP.p.1.clumped.ranges | gawk -F. '{ print $1,$2,$4 }' > REP.p.1.interval

# Run INRICH using all SNPs on msigdb 3.0 Pathways on REP set:
rm p1.paths
for rep1 in `cat top.1.paths.pre`
do
        grep -w $rep1 msigbd3.0.set >> p1.paths
done

nohup ./inrich \
-a REP.p.1.interval \
-m GIANT.map \
-g entrez.gene.map \
-t p1.paths \
-w 20 \
-r 100000 \
-q 5000 \
-i 20 \
-j 200 \
-p 1 \
-d .1 \
-o REP.p.1 &

###########################################################
# Examine replication power at different thresholds       #
###########################################################
# subtract top 53 lines and keep 536 following pathway lines

tail -n 1576 GIANT.p.001.out.inrich | head -n 536  > GIANT.p.001.pa
tail -n 1876 REP.p.001.out.inrich | head -n 536  > REP.p.001.pa

tail -n 7573 GIANT.p.01.out.inrich | head -n 514  > GIANT.p.01.pa
tail -n 6449 REP.p.01.out.inrich | head -n 536  > REP.p.01.pa

tail -n 23500 GIANT.p.1.out.inrich | head -n 536  > GIANT.p.1.pa
tail -n 24021 REP.p.1.out.inrich | head -n 536  > REP.p.1.pa

for i in `cat g`
do
        grep -w $i INRICH.p005.genes.csv >> gene.temp
done