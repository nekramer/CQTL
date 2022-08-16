#!/usr/bin/env Rscript

library(data.table)
args = commandArgs(trailingOnly = TRUE)

# args[1]...args[n]: alleleCount files from each sample


variants <- lapply(args[1:length(args)], fread, select = "variantID") |>
                rbindlist() |>
                unique() |>
                setorder()

fwrite(variants, file = "output/AI/variants.csv", quote = FALSE)