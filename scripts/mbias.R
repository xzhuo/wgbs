mbias_plot <- function(txt, pdf) {
    library("ggplot2")
    library("patchwork")
    df <- read.table(txt,header=F)
    colnames(df)<-c("sample","context","read","position","methylated_total","unmethylated_total","methylated_perc","coverage")
    plot <- ggplot(df, aes(x=position)) + geom_line(y=methylated_perc, color=context) + geom_line(y=methylated_total, color=context) + theme_bw()
    ggsave(pdf, plot)
}

mbias_plot(snakemake@input[[1]], snakemake@output[[1]])