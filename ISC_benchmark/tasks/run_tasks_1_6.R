#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
cmd_args <- c(file.path("ISC_benchmark", "run_ISC_benchmark.R"), args)
status <- system2("Rscript", args = cmd_args)
quit(status = status)
