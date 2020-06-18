library(readr)
library(dplyr)

cond <- "luma"
cnv_matrix <- read_tsv(paste0("data/", cond, "-cnvs.tsv"))
genes <- unlist(cnv_matrix[, 1])
cnv_matrix <- as.matrix(cnv_matrix[,-1])
rownames(cnv_matrix) <- genes

row_nas <- rowSums(is.na(cnv_matrix))
table(row_nas)

# row_nas
#   1   2   3   4   5   6   7   8   9  15  16  19  49 134 150 210 
# 514 103  35  49   4  43   1  13   6   1   1   1  41  10   1 392 
# Remove rows with equal or more than 134 NAs. Drop 403 rows
genes <- setdiff(genes,
                 names(which(row_nas > 133)))
cnv_matrix <- cnv_matrix[genes, ]

cnv_stats <- apply(cnv_matrix, 1, function(rr){
  rr <- rr[!is.na(rr)]
  return(list(no_na = length(rr),
              mean = mean(rr), 
              median = median(rr),
              sd = sd(rr),
              lowerq = quantile(rr, probs = 0.25),
              upperq = quantile(rr, probs = 0.75),
              min = min(rr), 
              max = max(rr)
              ))
})
cnv_stats = bind_rows(cnv_stats)
cnv_stats <- cnv_stats %>% mutate(ensemblID = genes) %>% select(ensemblID, everything())
cnv_stats <- cnv_stats %>% arrange(desc(sd))

write_tsv(cnv_stats, paste0("./data/", cond, "-genestats.tsv"))
