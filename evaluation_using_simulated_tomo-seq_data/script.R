library(tomoseqr)
library(SPsimSeq)
library(tidyverse)
rm(list=ls())

source("data/functions.R")

file_name = "qnap/20240810_five_genes/2000genes/result10.rda"
number_of_background <- 1995

# 1. background gene を voxel ごとにつくる

simulatedSeq_genes <- SPsimSeq(
  n.sim = 1,
  s.data = zhang.data.sub$counts,
  group = zhang.data.sub$MYCN.status,
  n.genes = number_of_background + 1,
  batch.config = 1,
  group.config = c(0.5, 0.5),
  tot.samples = 9282,
  pDE = 0.1,
  lfc.thrld = 0.5,
  result.format = "list",
  return.details = TRUE
)

simulated_count <- simulatedSeq_genes$sim.data.list[[1]]$counts %>% t()

matrix_background_genes <- matrix(0, nrow = 50^3, ncol = number_of_background)



exp1 <- generateExpressionPattern(50, 50, 50, func1, 100, 0)
exp2 <- generateExpressionPattern(50, 50, 50, func1, 100) +
  generateExpressionPattern(50, 50, 50, maskFunc, 10) +
  generateExpressionPattern(50, 50, 50, func1, -10)
exp3 <- generateExpressionPattern(50, 50, 50, func3, 100)
exp4 <- generateExpressionPattern(50, 50, 50, func4, 100)
exp5 <- generateExpressionPattern(50, 50, 50, func5, 100)
mask <- generateExpressionPattern(50, 50, 50, maskFunc, 1)
mask_df <- tomoseqr:::matrixToDataFrame(mask)

exp1_df <- tomoseqr:::matrixToDataFrame(exp1) %>% as_tibble() %>% rename(gene1 = value)
exp2_df <- tomoseqr:::matrixToDataFrame(exp2) %>% as_tibble() %>% rename(gene2 = value)
exp3_df <- tomoseqr:::matrixToDataFrame(exp3) %>% as_tibble() %>% rename(gene3 = value)
exp4_df <- tomoseqr:::matrixToDataFrame(exp4) %>% as_tibble() %>% rename(gene4 = value)
exp5_df <- tomoseqr:::matrixToDataFrame(exp5) %>% as_tibble() %>% rename(gene5 = value)
matrix_background_genes[mask_df$value == 1,] %>% dim()
matrix_background_genes[mask_df$value == 1,] <- simulated_count
tibble_background_genes <- as_tibble(matrix_background_genes, .name_repair = "unique")
colnames(tibble_background_genes) <- paste0("gene", 6:(number_of_background + 5))

binded_each_genes <- exp1_df %>%
  full_join(exp2_df, by = join_by(x,y,z)) %>%
  full_join(exp3_df, by = join_by(x,y,z)) %>%
  full_join(exp4_df, by = join_by(x,y,z)) %>%
  full_join(exp5_df, by = join_by(x,y,z)) %>%
  bind_cols(tibble_background_genes)

sum_of_each_voxel <- binded_each_genes %>% select(-x, -y, -z) %>% apply(1, sum)
ratio <- 1e04 / sum_of_each_voxel
ratio[is.infinite(ratio)] <- 0

for (i in 4:length(binded_each_genes[1,])) {
    binded_each_genes[,i] <- binded_each_genes[,i] * ratio
}

indexList <- expand.grid(1:50, 1:50, 1:50)
exp1_aligned <- exp1
for (i in seq_along(indexList[,1])) {
  x <- indexList[i, 1]
  y <- indexList[i, 2]
  z <- indexList[i, 3]
  exp1_aligned[x, y, z] <- binded_each_genes$gene1[i]
}
exp2_aligned <- exp2
for (i in seq_along(indexList[,1])) {
  x <- indexList[i, 1]
  y <- indexList[i, 2]
  z <- indexList[i, 3]
  exp2_aligned[x, y, z] <- binded_each_genes$gene2[i]
}
exp3_aligned <- exp3
for (i in seq_along(indexList[,1])) {
  x <- indexList[i, 1]
  y <- indexList[i, 2]
  z <- indexList[i, 3]
  exp3_aligned[x, y, z] <- binded_each_genes$gene3[i]
}
exp4_aligned <- exp4
for (i in seq_along(indexList[,1])) {
  x <- indexList[i, 1]
  y <- indexList[i, 2]
  z <- indexList[i, 3]
  exp4_aligned[x, y, z] <- binded_each_genes$gene4[i]
}
exp5_aligned <- exp5
for (i in seq_along(indexList[,1])) {
  x <- indexList[i, 1]
  y <- indexList[i, 2]
  z <- indexList[i, 3]
  exp5_aligned[x, y, z] <- binded_each_genes$gene5[i]
}


tomo_X <- binded_each_genes %>% group_by(x) %>% summarize_all(.funs = sum) %>% select(-x,-y,-z) %>% t() %>% as_tibble(.name_repair = "minimal")
tomo_Y <- binded_each_genes %>% group_by(y) %>% summarize_all(.funs = sum) %>% select(-x,-y,-z) %>% t() %>% as_tibble(.name_repair = "minimal")
tomo_Z <- binded_each_genes %>% group_by(z) %>% summarize_all(.funs = sum) %>% select(-x,-y,-z) %>% t() %>% as_tibble(.name_repair = "minimal")
gene_id_list <- colnames(binded_each_genes)[4:length(binded_each_genes)]
section_list <- paste0("section", 1:50)
colnames(tomo_X) <- section_list
tomo_X <- bind_cols(tibble(gene_id =gene_id_list), tomo_X)
colnames(tomo_Y) <- section_list
tomo_Y <- bind_cols(tibble(gene_id =gene_id_list), tomo_Y)
colnames(tomo_Z) <- section_list
tomo_Z <- bind_cols(tibble(gene_id =gene_id_list), tomo_Z)

withNormalization <- estimate3dExpressions(
  tomo_X,
  tomo_Y,
  tomo_Z,
  mask,
  query = c("gene1", "gene2", "gene3", "gene4", "gene5"),
  normalize = TRUE
)

withoutNormalization <- estimate3dExpressions(
  tomo_X,
  tomo_Y,
  tomo_Z,
  mask,
  query = c("gene1", "gene2", "gene3", "gene4", "gene5"),
  normalize = FALSE
)

groundtruth <- withNormalization
groundtruth$results$gene1$reconst <- exp1_aligned
groundtruth$results$gene2$reconst <- exp2_aligned
groundtruth$results$gene3$reconst <- exp3_aligned
groundtruth$results$gene4$reconst <- exp4_aligned
groundtruth$results$gene5$reconst <- exp5_aligned

save(binded_each_genes, tomo_X, tomo_Y, tomo_Z, mask, withNormalization, withoutNormalization, groundtruth, file = file_name)
