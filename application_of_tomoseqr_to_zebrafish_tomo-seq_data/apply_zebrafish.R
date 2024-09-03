library(tidyverse)
library(tomoseqr)

# Mask generator
generateMask <- function(
    xLen,
    yLen,
    zLen,
    func
) {
  mask <- array(0, dim = c(xLen, yLen, zLen))
  indexList <- expand.grid(1:xLen, 1:yLen, 1:zLen)
  for (i in seq_along(indexList[, 1])) {
    x <- indexList[i, 1]
    y <- indexList[i, 2]
    z <- indexList[i, 3]
    if (func(x, y, z)) {
      mask[x, y, z] <- 1
    }
  }
  return(mask)
}

# Definition of the shape of the mask
maskFunction <- function(x, y, z) {
  return(
    (x - 33)^2 + (y - (49 / 2))^2 + (z - (56 / 2))^2 <= (44 / 2)^2 &
      (x - 33)^2 + (y - (49 / 2))^2 + (z - (56 / 2))^2 >= (17)^2 &
      x <= 37
  )
}


# Data loading
zebra_AV <- read_csv("data/zebra_AV.csv")
zebra_VD <- read_csv("data/zebra_VD.csv")
zebra_LR <- read_csv("data/zebra_LR.csv")
zebra_gene_list <- read_csv("data/zebra_gene_list.csv") %>% pull(geneID)
mask <- generateMask(50, 49, 56, func = maskFunction)

# Reconstruction
zebra_reconst_result <- estimate3dExpressions(
  zebra_AV,
  zebra_VD,
  zebra_LR,
  mask,
  query = zebra_gene_list
  )
