MakeColorKey <- function(break_seq,heat_map_colors,save_directory)
#Makes a color key for heatmaps because the default one is a little off (the first segment is too thin and the last too thick)
{
n_colors <- length(break_seq)-1
interval <- break_seq[2]-break_seq[1]
x <- seq(break_seq[1]+0.5*interval,tail(break_seq,n=1)-0.5*interval,interval)
pdf(paste(save_directory,'color_key.pdf'),width=5,height=3.5)
par(yaxs="i") #removes gap between x-axis and bars
hist(x, col=heat_map_colors,breaks=break_seq, border=NA,xlab=NA,ylab=NA,main=NA)
axis(side = 1, lwd=1)
dev.off()
}