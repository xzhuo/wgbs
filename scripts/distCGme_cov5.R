#!/opt/htcf/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/r-3.5.1-l2myjcnyryfxpteoeu4zbqcigtpm6mtz/bin/Rscript
###!/opt/apps/r/3.3.3/bin/Rscript

# Author: Hyung Joo Lee

### library, working directory
library(ggplot2)

### input
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
wfile <- args[2]

### output
base <- gsub("^.*/", "", file)
base <- gsub(".dss.txt.gz", "", base)
#wfile <- paste0(base, "_CGme_density.txt")

## Load file
cpme <- read.delim(file = file, as.is = TRUE, header=T)
cpme <- cpme[ which(cpme[,1]!="chrM"), ]
cpme <- cpme[, 3:4]

## Coverage at least 5
cpme <- cpme[ which(cpme[,1] >= 5), ]

## Calculate % CpGs
cpme[,1] <- (cpme[,2] / cpme[,1]) *100

### density graph
g <- ggplot(cpme, aes(cpme[,1])) + geom_density(adjust=5, n=101)

### plot coordiantes
g <- ggplot_build(g)
density <- c(base, g$data[[1]][,1])
write(density, wfile, sep ="\t")

q()
