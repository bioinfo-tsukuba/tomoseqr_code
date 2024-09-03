library(tidyverse)
library(tomoseqr)
# Data loading
pla_AP <- read_csv("data/pla_AP.csv")
pla_LR <- read_csv("data/pla_LR.csv")
pla_DV <- read_csv("data/pla_DV.csv")
pla_gene_list <- read_csv("data/pla_gene_list.csv") %>% pull(geneID)
pla_gene_table <- pla_gene_list %>% as_tibble_col() %>% mutate("geneName" = c("piwiA", "opsin", "Djf-1", "DjNp19"))
load("data/pla_mask.rda")

# Reconstruction
pla_reconst_result <- estimate3dExpressions(
  pla_AP,
  pla_LR,
  pla_DV,
  mask,
  query = pla_gene_list,
  normalize = FALSE
)

reconstruct_top_5000 <- estimate3dExpressions(
  pla_AP,
  pla_LR,
  pla_DV,
  mask,
  query = top_5000_genes,
  normalize = FALSE
)
