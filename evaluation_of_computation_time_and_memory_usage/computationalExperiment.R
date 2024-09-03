library(tomoseqr)
library(readr)

AV <- read_csv("data/zebra_AV.csv", show_col_types = FALSE)
VD <- read_csv("data/zebra_VD.csv", show_col_types = FALSE)
LR <- read_csv("data/zebra_LR.csv", show_col_types = FALSE)
load("data/zebra_mask.rda")

geneList <- extractGeneList(AV, VD, LR)
system.time({result <- estimate3dExpressions(AV, VD, LR, mask = mask, query = geneList[1:100])})
