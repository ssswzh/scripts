library(data.table)
library(reshape2)

usage <- "
  Usage:
    Rscript --vanilla draw_wgs_depth_distribution.r samtools.depth.tsv fasta.bed output.name
    
  samtools.depth.tsv:
    samtools depth sample.bam > samtools.depth.tsv
  
  fasta.bed (TAB-delimited):
    chr1	1	141
    chr2	1	304
  
"

args <- commandArgs(TRUE)
if (length(args)==0) {
  stop(usage, call.=FALSE)
} else if (args[1]=='-h'||args[1]=='--help') {
  stop(usage, call.=FALSE)
}

depth <- args[1]
chrom <- args[2]
output <- args[3]

windows = 100000

df = fread(depth,header=F,sep="\t")
df = as.data.frame(df)
df[,2] = ceiling(df[,2]/windows)
m = melt(df,id=c("V1","V2"))
df=dcast(m,V1+V2~variable,mean)
#df$V4 = df$V3/max(df$V3)
df$V4 = sqrt(df$V3)/max(sqrt(df$V3))
maxdepth = ceiling(max(df[,3]))
#if(maxdepth>1000) region = c(50,maxdepth/10:1) else region = maxdepth/10:1
region = rev((sqrt(maxdepth)/1:12)^2) ; if(region[1]>200) region=c(50,region)
colindex = c()
for(j in 1:length(df[,3])){
	for(i in 1:length(region)){
	if(i==1){
		if(df[j,3]<=region[1])
		colindex = c(colindex,i)
	}else{
		if(df[j,3]>region[i-1]&df[j,3]<=region[i])
		colindex = c(colindex,i)
	}
	}
}
col = c("darkgreen", "yellow", "red")
col = colorRampPalette(col)(length(region))
color = col[colindex]
df$V5 = color
print(table(color))
## chrom
t = read.table(chrom,header=F,sep="\t")
t[,3] = ceiling(t[,3]/windows)
nchr = nrow(t)

pdf(output)
chr_width = 1 ; d = 0.2
y1 = 0:(nchr-1)*(-(chr_width+d))
y2 = y1 -1
par(mgp=c(3,0.5,0))
plot(1,xlim=c(0,max(t[,3])),ylim=c(min(y2)-d,0),type="n",axes=F,main="Reads Distribution(Kb)",xlab="",ylab="")
rect(0,y2,t[,3],y1,col="grey70",border="grey70")
par(xpd=T)
text(0,(y1+y2)/2,t[,1],pos=2,cex=0.8)
for(i in 1:nrow(t)){
dfi = df[df[,1]==t[i,1],,drop=F]
if(nrow(dfi)>0)segments(dfi[,2],y2[i],dfi[,2],y2[i]+dfi[,4],col=dfi[,5])
}
axis(3)
legend("bottomright",legend=ceiling(region),pt.cex=3,pch=15,col=col,bty="n")
dev.off()
