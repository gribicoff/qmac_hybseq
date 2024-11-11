library(tidyverse)

# read in dataframes of number of SNPs per gene, functional class annotation for each gene

nsnp.df <- read.table("numsnpspergene.txt", header = TRUE) %>%
    dplyr::rename(genefull = gene) %>%
    dplyr::mutate(gene = gsub("_.*", "", genefull))
class.df <- read.table("allgenes_ids.txt", header = TRUE)

# join class with number of SNPs per gene

gene.df <- left_join(nsnp.df, class.df, by = "gene") %>%
    dplyr::select(-gene) %>%
    dplyr::filter(n_snps >= 1)

# write out list of genes per functional class with at least one SNP

classes <- c("other", "bud", "drought")
sapply(classes, function(x) {
    gene.df %>%
        filter(class == x) %>%
        dplyr::select(genefull) %>%
        write.table(paste0(x, "genes_idsforsub.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
})
