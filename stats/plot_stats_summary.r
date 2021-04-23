#!/usr/bin/env Rscript

# Rscript plot_stats_summary.r file_name column column_name

args = commandArgs(trailingOnly=TRUE)
if(length(args)!=3){
    stop("Usage: Rscript plot_stats_summary.r file_name column column_name")
}

library(fBasics)
data <- read.table(args[1], header=F)
#sub_data <- subset(data, V1 < 30000)
#sub_data <- data
colu <- as.integer(args[2])
if(args[3]=="Read_identity"||args[3]=="Read_Quality"){data[,colu] <- data[,colu]*100}
if(max(data[,colu])-min(data[,colu])>500){breaks <- 500} else if(max(data[,colu])-min(data[,colu])>50){breaks <- seq(min(data[,colu])-0.5,max(data[,colu])+1.5,1)} else{breaks <- seq(min(data[,colu])-0.5,max(data[,colu])+0.5,0.1)}

png(filename=paste0(args[1],".",args[3],".png"), width=1280, height=960, units="px")
#h <- hist(data[,colu], xlab=args[3], main=paste0(basename(args[1]),".",args[3]), breaks=breaks, cex.axis=1.5, cex.lab=2, cex.main=2)
h <- hist(data[,colu], xlab=args[3], main=args[3], breaks=breaks, cex.axis=1.5, cex.lab=2, cex.main=2)
x_line <- h$mid[which(h$counts == max(h$counts))]
abline(v=x_line, col="red", lwd=1.5, lty=2)
text(x_line, max(h$counts), paste0("Max:",x_line), col="red")
box()
dev.off()

out <- data.frame(mid = h$mid, counts = h$counts, density = h$density)
write.table(out, file = paste0(args[1],".hist_",args[3],".txt"), row.names = FALSE, sep = "\t", quote = FALSE)

write.table(basicStats(data[,colu]), file = paste0(args[1],".summary_",args[3]), sep = "\t", quote = FALSE)
