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
base <- gsub(".tmp.insert.txt.gz", "", base)


## Load file
insert <- read.delim(file = file, as.is = TRUE, header=F)


### density graph
g <- ggplot(insert, aes(insert[,4])) + geom_density(n=126) + xlim(c(0,500))


### plot coordiantes
g <- ggplot_build(g)
density <- c( base, g$data[[1]][,1] )
write(density, wfile, sep ="\t")

q()
