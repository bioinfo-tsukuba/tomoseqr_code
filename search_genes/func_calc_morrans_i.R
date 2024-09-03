library(spdep)
library(tidyverse)
library(tomoseqr)

calc_morrans_i <- function(tomoObj, geneID) {
  if(geneID %in% names(tomoObj[["results"]]) == FALSE){
    return(c(NA, NA))
  }
  if(var(tomoObj[["results"]][[geneID]][["reconst"]]) == 0) {
    return(c(NA, NA))
  }
  dims <- dim(tomoObj[["results"]][[geneID]][["reconst"]])
  coords <- tomoseqr:::matrixToDataFrame(tomoObj[["mask"]]) %>%
    filter(value != 0) %>%
    select(-value)
  nb <- knn2nb(knearneigh(coords, k = 26))
  listw <- nb2listw(nb, style = "W")
  filterd_gene_data <- coords %>% left_join(toDataFrame(tomoObj, geneID))
  gene_data_flat <- filterd_gene_data %>% pull(value)
  moran_result <- moran.test(gene_data_flat, listw)
  
  return(c(as.numeric(moran_result[["estimate"]][1]), moran_result[["p.value"]]))
}
