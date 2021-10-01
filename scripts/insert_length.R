insert_len_plot <- function(txt, pdf) {
	library("ggplot2")
	df <- read.table(txt,header=T)
	plot <- ggplot(df, aes(x=position, y=density,color=sample)) + geom_line(alpha=0.5) + xlab("Insertion length (bp)") +ylab("Density") + theme_bw()
	ggsave(pdf, plot)
}

insert_len_plot(snakemake@input[[1]], snakemake@output[[1]])