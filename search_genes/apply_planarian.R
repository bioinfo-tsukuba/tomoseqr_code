rm(list = ls())
library(tidyverse)
library(tomoseqr)
# Data loading
pla_AP <- read_csv("data/pla_AP.csv")
pla_LR <- read_csv("data/pla_LR.csv")
pla_DV <- read_csv("data/pla_DV.csv")
exp_descending <- pla_AP %>% select(-geneID) %>% apply(1, sum) %>% as_tibble_col() %>% bind_cols(pla_AP %>% select(geneID)) %>% arrange(-value)
gene_IDs <- exp_descending[20001:25566,] %>% pull(geneID)
gene_IDs <- c(gene_IDs, "DjGI008464_001")
load("data/pla_mask.rda")

# Reconstruction
reconstruct_result <- estimate3dExpressions(
  pla_AP,
  pla_LR,
  pla_DV,
  mask,
  query = gene_IDs,
  normalize = FALSE
)

save(reconstruct_result, file = "qnap/pla_5000/result_20001_25566.rda")
