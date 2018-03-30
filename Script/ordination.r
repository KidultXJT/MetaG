Args <- commandArgs()

# Necessary
intable      <- Args[6]
outDir       <- Args[7]
prefix       <- Args[8]

library(Kitult)
require(pheatmap)
require(ggplot2)
require(methods)

