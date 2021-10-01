preseq_plots <- function(txt, pdf) {
    library("ggplot2")
    df <- read.table(txt,header=T)
    plot <- ggplot(df, aes(x=TOTAL_READS, y=EXPECTED_DISTINCT,color=SAMPLE)) + geom_abline(intercept = 0, slope = 1)+ geom_line() + xlim(c(0, 5.0e+8)) + ylim(c(0, 5.0e+8)) + xlab("Total mapped reads") +ylab("Expected distinct reads") + theme_bw()
    ggsave(pdf, plot)
}

preseq_plots(snakemake@input[[1]], snakemake@output[[1]])