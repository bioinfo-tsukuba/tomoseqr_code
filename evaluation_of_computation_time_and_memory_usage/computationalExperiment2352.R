library(tomoseqr)
library(readr)
library(dplyr)
library(tibble)

df_to_array <- function(df) {
  x <- df[["x"]]
  y <- df[["y"]]
  z <- df[["z"]]
  ret_array <- array(0, dim = c(max(x), max(y), max(z)))
  mask_val <- df[["value"]]
  for (i in seq_along(x)) {
      ret_array[x[i], y[i], z[i]] <- mask_val[i]
  }
  return(ret_array)
}

AV <- read_csv("data/zebra_AV.csv", show_col_types = FALSE)
VD <- read_csv("data/zebra_VD.csv", show_col_types = FALSE)
LR <- read_csv("data/zebra_LR.csv", show_col_types = FALSE)
load("data/zebra_mask.rda")

AVmax <- tibble(GeneID=AV$GeneID, max=apply(AV[, -1], 1, FUN=max)) %>%
  filter(max > 50)
geneList <- intersect(extractGeneList(AV, VD, LR), AVmax$GeneID)
system.time({result <- estimate3dExpressions(AV, VD, LR, mask = mask, query = geneList)})
