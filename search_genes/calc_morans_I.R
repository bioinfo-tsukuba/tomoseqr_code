rm(list=ls())
library(tidyverse)
library(tomoseqr)
source("data/func_calc_morrans_i.R")
load("data/Djb5_2mouse_entrez")
go_anno <- read_tsv("data/QuickGO-annotations-1724821892765-20240828.tsv")

pla_id <- Dj2mouse %>% filter(GeneName %in% go_anno$SYMBOL) %>% pull(planarian)

pla_AP <- read_csv("data/pla_AP.csv")
pla_LR <- read_csv("data/pla_LR.csv")
pla_DV <- read_csv("data/pla_DV.csv")
load("data/pla_mask.rda")

# Reconstruction
reconstruct_result <- estimate3dExpressions(
  pla_AP,
  pla_LR,
  pla_DV,
  mask,
  query = pla_id,
  normalize = FALSE
)

vec_calc_moran <- Vectorize(calc_morrans_i, vectorize.args = "geneID")
morans_i <- vec_calc_moran(tomoObj = reconstruct_result, geneID = pla_id)
morans_df <- tibble(geneID=colnames(morans_i), morans_index=morans_i[1,], pvalue=morans_i[2,])
write_csv(morans_df, file = "data/morans_i_df.csv")

morans_df %>% ggplot(aes(x=morans_index)) +
  geom_histogram() +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank())

morans_top_5 <- morans_df %>%
  arrange(-morans_index) %>%
  filter(row_number() <= 5) %>%
  pull(geneID)

reconstruct_result_top5 <- estimate3dExpressions(
  pla_AP,
  pla_LR,
  pla_DV,
  mask,
  query = morans_top_5,
  normalize = FALSE
)
