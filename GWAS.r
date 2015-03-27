# GWAS:

system("plink --bfile CARDIA.clean.FINAL --pheno CARDIA.clean.FINAL.phe --linear --out CARDIA.clean.FINAL")
system("plink --bfile ARIC.clean.FINAL --pheno ARIC.clean.FINAL.phe --linear --out ARIC.clean.FINAL")
system("plink --bfile MESA.clean.FINAL --pheno MESA.clean.FINAL.phe --linear --out MESA.clean.FINAL")

system("plink --bfile MERGE.clean.FINAL --pheno MERGE.clean.FINAL.phe --linear --out MERGE.clean.FINAL")
