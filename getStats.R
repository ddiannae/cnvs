library(readr)
library(dplyr)

cond <- "basal"
cnv_matrix <- read_tsv(paste0("data/", cond, "-cnvs.tsv"))
genes <- unlist(cnv_matrix[, 1])
cnv_matrix <- as.matrix(cnv_matrix[,-1])
rownames(cnv_matrix) <- genes

row_nas <- rowSums(is.na(cnv_matrix))
table(row_nas)

# row_nas luma
#   1   2   3   4   5   6   7   8   9  15  16  19  49 134 150 210 
# 514 103  35  49   4  43   1  13   6   1   1   1  41  10   1 392 
# Remove rows with equal or more than 134 NAs. Drop 403 rows

# row_nas lumb
# 0         1     2     3     4     5     6     7     8     9    11    22    66   133   159   189 
# 58697   953   169    52    67     2    41     4     1     2    13     1    41    10     1   392
# Remove rows with equal or more than 133 NAs. Drop 403 rows

# row_nas her2
#     0     1     2     3     4     5    54    72    77   101 
# 58912   856   187    15    31     1    41     1    10   392
# Remove rows with equal or more than 54 NAs. Drop 444 rows

# row_nas basal
#     0     1     2     3     4     5     6     7     8    11    12    21   131   157   160   215 
# 58483  1188   177    76    57     2    12     1     1     2     2     1    41     1    10   392
# Remove rows with equal or more than 131 NAs. Drop 444 rows

genes <- setdiff(genes,
                 names(which(row_nas > 130)))
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
cnv_stats <- cnv_stats %>% mutate(ensembl_id = genes) %>% 
  select(ensembl_id, everything()) %>% 
  arrange(desc(sd))

write_tsv(cnv_stats, paste0("./data/", cond, "-genestats.tsv"))
