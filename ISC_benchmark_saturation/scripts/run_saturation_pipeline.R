#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
})

targets::tar_make()
