#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("\nUsage:
  Rscript --vanilla plot_fragment_in_genome_region.R iris.txt out.txt
  ", call.=FALSE)
}

library("optparse")
 
option_list = list(
  make_option("--file", type="character", default=NULL, help="dataset file name", metavar="character"),
  make_option("--prefix", type="character", default="out.txt", help="output file prefix [default= ${filename}.pdf]", metavar="character"),
  make_option("--start", type="character", default="Alignment_start", help="column name of fragment start [default= %default]", metavar="character"),
  make_option("--end", type="character", default="Alignment_end", help="column name of fragment end [default= %default]", metavar="character"),
  make_option("--sample", type="character", default="ReadID", help="column name of sample name [default= %default]", metavar="character"),
  make_option("--pos", type="character", default=NULL, help="positions to draw vertical lines in plot", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


library(ggplot2)
library(ggalt)
theme_set(theme_classic())
ma <- read.table(opt$file, header=T, sep="\t", row.names=1)
ma$Sample <- factor(rownames(ma), levels=as.character(rownames(ma)))
gg <- ggplot(ma, aes(x=opt$start, xend=opt$end, y=opt$sample)) + 
  geom_dumbbell(colour_x="#0e668b", color="#a3c4dc", size=1.5, colour_xend="black") + 
  #scale_x_continuous(label=percent) + 
  labs(x="Coordinate", y="Contigs", title="PSDL-2 Alpha Haplotype2 in Hg19 Position") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.background=element_rect(fill="#f7f7f7"),
        panel.background=element_rect(fill="#f7f7f7"),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(),
        axis.ticks=element_blank(),
        legend.position="top",
        panel.border=element_blank()) 
for(i in opt$pos){
  gg <- gg + geom_vline(xintercept = i,color="red")
}

pdf(opt$prefix+".pdf")
plot(gg)
dev.off()
