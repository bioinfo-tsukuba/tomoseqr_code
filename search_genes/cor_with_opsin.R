# By changing RECONSTRUCTION_RESULT_PATH and executing it repeatedly, compile the analysis results of all files into a single CSV file.

rm(list = ls())
library(tidyverse)
library(tomoseqr)

RECONSTRUCTION_RESULT_PATH <- "qnap/pla_5000/result_20001_25566.rda"

opsin_id <- "DjGI008464_001"
load(RECONSTRUCTION_RESULT_PATH)
other_result <- read_csv("data/cor_result.csv")
cor_result <- correlationWithSpecificGene(reconstruct_result, opsin_id)
cor_result <- cor_result %>% bind_rows(other_result)
write_csv(cor_result, file = "data/cor_result.csv")
