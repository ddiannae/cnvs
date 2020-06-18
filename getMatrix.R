library(readr)
library(dplyr)
library(stringr)

cond <- "luma"
files_to_read <- list.files(path = paste0("./raw/", cond), pattern = "\\.tsv$", full.names = T)

all_files <- lapply(files_to_read, function(file) {
  cnvs <- read_tsv(file)
  cnvs <- cnvs %>% mutate(gene_id = as.character(lapply(str_split(gene_id, "\\."), "[[", 1)))
  return(cnvs %>% select(gene_id, copy_number))
})
  
size <- unique(do.call(rbind,lapply(all_files, dim)))
stopifnot(nrow(size)==1)

genes <- do.call(cbind, lapply(all_files, function(x) select(x, "gene_id")))
genes <- t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))

targets <- data.frame(id = paste(cond, 1:length(files_to_read), sep = "_"), 
                      file = str_remove(str_remove(lapply(str_split(files_to_read, "/"), "[[", 4), 
                                                   "TCGA-BRCA."), ".gene_level_copy_number.tsv"))


matrix <- bind_cols(lapply(all_files, function(x) select(x, "copy_number")))
colnames(matrix) <- targets$id
matrix <- matrix %>% mutate(ensembl_id = genes[,1]) %>% 
  select(ensembl_id, everything())

write_tsv(matrix, paste0("./data/", cond, "-cnvs.tsv"))
write_tsv(targets, paste0("./data/", cond, "-files.tsv"))
