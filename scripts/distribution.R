distribution <- function(txt, pdf) {
	library("ggplot2")
	df <- read.table(txt,header=T)
	plot <- ggplot(df, aes(x=perc,y=count,color=sample)) + geom_smooth(se = FALSE, size =0.5, alpha=0.2) + ylab("density")+ theme_bw()
	ggsave(pdf, plot)
}

distribution(snakemake@input[[1]], snakemake@output[[1]])
