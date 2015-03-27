# Function Inputs:
# 1. name of dataset #
# 2. up/down stream additional base positions beyond gene to include as genic
# 3. number of cores to use
# 4. name of file containing list of unique autosomal genes with header (NAME, START, STOP, CHR) ; NOTE: ALL CHARACTERS IN GENE NAMES MUST BE CAPITALIZED
#
# Function Outputs:
# 1. gene set file (.set extension)
# 2. .bim, .bed, .fam files for each gene

gene_snp_sets <- function(dataset,KB,CORES,clean.genelist){ # dataset and clean.list must be in quotes

  total.snps <- read.table(paste('',dataset,'.bim',sep="")) # read in all SNPs
  names(total.snps) <- c('CHR','SNP','CM','BP','A1','A2')
  
  gene.locations <- read.table(clean.genelist,header=TRUE) # read in all genes

  snp.info.data <- as.data.frame(cbind(as.character(total.snps$SNP),total.snps$CHR,total.snps$BP))
  gene.info.data <- as.data.frame(cbind(as.character(gene.locations$NAME),gene.locations$CHR,gene.locations$START,gene.locations$STOP))
  

    load(paste('',dataset,'.gene.set.Rdata',sep="")) # read in gene set information created using 'generate_gene_set()' function
  path.genes <- list()
  for(i in 1:512){
    path.genes[i] <- as.list(gene.sets.final[[i]])
    print(i)
  }

  path.lengths <- unlist(lapply(path.genes,length))

  
