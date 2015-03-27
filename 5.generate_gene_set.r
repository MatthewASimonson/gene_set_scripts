# Function Inputs:
# 1. name of .csv file that contains pathway information (1 pathway per line; expecting msigdb format)
# 2. name of dataset
# 3. file with list of genes being examined; format of header ("NAME","CHR","START","STOP"), 1 gene per row
# 4. mininum number of genes allowed per set (20 genes)
# 5. maximum number of genes allowed per set (200 genes) Wang et al. (2007)
#
# Function Outputs:
# 1. .Rdata file containing object with gene set information (GSEA.set) in proper format for gseaGene() function; genes not included in input file will be removed and not counted


generate_gene_set <- function(csv_name,dataset,gene.info,min_gene,max_gene){

# Read in pathway information (msigdb format)
#############################################  
pathways <- read.table(csv_name,header=FALSE,sep=",",fill=TRUE) #

set_assign <- function(row){ # Function that reads genes listed in pathways data.frame of gene names (1 path per row)
     temp <- as.matrix(row[3:length(row)])
     keep.index <- which(temp!="")
     GENES <- temp[keep.index]
     genes.in.path <- as.data.frame(GENES)
     return(genes.in.path)
   } # END set_assign

gene.sets <- apply(pathways, 1, set_assign) # apply set_assign function to bin.genes list of genes data.frame
names(gene.sets) <- pathways$V1[1:length(gene.sets)]

gene.list <- read.table(gene.info,header=TRUE) # read in list of genes; "clean.genelist"

check_path <- function(g.set){ # Function that scans all pathways and drops any genes not included in gene.list
     g.set <- as.data.frame(g.set)
     g.matches <- as.matrix(g.set) %in% gene.list$NAME # which genes in pathway are included in gene list
     GENES <- g.set$GENES[g.matches]
     out.set <- as.data.frame(GENES)
     return(out.set)
} # END check_path

check.sets <- lapply(gene.sets, check_path) # drop any genes listed in pathways not included in gene list (run time = 22 seconds approx)

remove_min_max <- function(check.sets,min_gene,max_gene){ # function removes any pathways that have too many or too few genes
     # return pathway list with NA's in pathways with too many or too few genes:
     min_max_paths <- function(na.sets,min_gene,max_gene){ 
       na.sets <- as.data.frame(na.sets)
       if ((nrow(na.sets)<min_gene) | (nrow(na.sets)>max_gene)){
         na.sets <- NA
       }
       return(na.sets)  
       } # END min_max_paths

  na.sets <- lapply(check.sets,min_max_paths,min_gene,max_gene) # sets with NA's in excluded indeces
  out.sets <- na.sets[-c(which(is.na(na.sets)))]
  return(out.sets)   
} # END remove_min_max

gene.sets.final <- remove_min_max(check.sets,min_gene,max_gene) # remove pathways that have too many or too few genes

save(gene.sets.final,file=paste('',dataset,'.gene.set.Rdata',sep="")) # save object to current directory
}  # END check_list

csv_name <- "canonical_path.csv"
dataset <- 'MERGE.clean.FINAL'
gene.info <- "clean.genelist.FINAL" 
min_gene <- 20
max_gene <- 200
generate_gene_set(csv_name,dataset,gene.info,min_gene,max_gene)
