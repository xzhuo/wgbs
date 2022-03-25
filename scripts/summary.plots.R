aggregate_plots <- function(alignment, conversion, reads_number, reads_perc, conversion_pdf) {
    library("ggplot2")
    library("dplyr")
    library("tidyr")
    df <- read.table(alignment,header=T)
    reads <- df %>% select (sample, total, after_trim, uniq, dedup) %>% gather(type,reads,total,after_trim, uniq, dedup)
    reads$type <- factor(reads$type, levels = c("total","after_trim","uniq","dedup"))
    reads <- ggplot(reads,aes(x=sample,y=reads,fill=type))+geom_bar(stat="identity",position=position_dodge())+coord_flip()+theme_bw()
    ggsave(reads_number, reads)
    perc <- df %>% select (sample, unique_perc, dedup_perc) %>% gather(type,percentage,unique_perc,dedup_perc)
    perc$percentage <- as.numeric(sub("%","",perc$percentage)) / 100
    perc <- ggplot(perc,aes(x=sample,y=percentage,fill=type))+geom_bar(stat="identity",position=position_dodge())+scale_y_continuous(labels = scales::percent_format())+coord_flip(ylim = c(0.5,1))+theme_bw()
    ggsave(reads_perc, perc)

    df <- read.table(conversion,header=T)
    if ("CH_conv_rate" %in% colnames(df)) {
        df <- df %>% select (sample, lambda_rate, CH_conv_rate) %>% gather(type,conversion_rate,lambda_rate,CH_conv_rate)
    } else {
        df <- df %>% select (sample, lambda_rate) %>% gather(type,conversion_rate,lambda_rate)
    }
    
    df$conversion_rate <- as.numeric(sub("%","",df$conversion_rate)) / 100
    gg <- ggplot(df,aes(x=sample,y=conversion_rate,fill=type))+geom_bar(stat="identity",position=position_dodge())+scale_y_continuous(labels = scales::percent_format())+coord_flip(ylim = c(0.98,1))+theme_bw()
    ggsave(conversion_pdf, gg)
}

aggregate_plots(snakemake@input[[1]], snakemake@input[[2]], snakemake@output[[1]], snakemake@output[[2]], snakemake@output[[3]])