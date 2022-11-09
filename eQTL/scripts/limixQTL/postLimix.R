#!/usr/bin/R
library(rhdf5)
args <- commandArgs(trailingOnly = TRUE)

# List the content of the file
dtlist <- h5ls(args[1])$name