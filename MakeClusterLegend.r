MakeClusterLegend <- function(n_clusters,cluster_colors,save_directory)
{


            cluster_names <- 1:n_clusters
            x=rep(0.25,length(cluster_names))
            pdf(paste(save_directory,'cluster_color_legend','.pdf',sep=''),width=5,height=5)
            barplot(x, col=cluster_colors,border=NA,xlab=NA,ylab=NA,ylim=c(0,4))
            legend(x=1,y=4,legend=cluster_names,col=cluster_colors,pch = c(15, 15, 15, 15),xpd=NA,bty='n', inset=c(1,1),pt.cex=2)
            dev.off()

}
