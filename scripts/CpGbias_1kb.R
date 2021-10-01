#!/opt/htcf/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/r-3.5.1-l2myjcnyryfxpteoeu4zbqcigtpm6mtz/bin/Rscript
##!/opt/apps/r/3.3.3/bin/Rscript
# Author: Hyung Joo Lee

### library, working directory
library(ggplot2)


### INPUT
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
name <- args[2]


### OUTPUT
wfile <- paste0(name, ".CpG_bias.boxplot.pdf")


## Load file
comp <- data.frame(CpG=numeric(), Coverage=numeric())
comp <- read.delim(file,
                   col.names = c("", "", "", "CpG", "Coverage"))[,c(4,5)]

print("Pearson, coverage")
cs<-cor.test(comp$CpG, comp$Coverage,method="pearson")
cs
cs$estimate

print("Spearman")
cs<-cor.test(comp$CpG, comp$Coverage,method="spearman")
cs
cs$estimate


## Plot the boxplot
pdf(wfile, 3,5)

g <- ggplot(comp, aes(as.factor(CpG), Coverage) ) 
g + geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(name = "CpG counts / 100 bp",breaks=c(1,5,10,15,20))+
  scale_y_continuous(name = "Read coverage")+
  coord_flip(ylim=c(0,100)) +
  theme_bw() +
  theme(text = element_text(size=16, color="black"),
        axis.text=element_text(size=14, color="black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.box.background=element_blank(),
        legend.key=element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
  )

dev.off()

q()

