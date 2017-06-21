MakeGroupLegend <- function(groups_corresponding,GroupColorMatrix,ColGroupsScheme,save_directory)
{
    if (!is.null(GroupColorMatrix))
    {
        n_legends <- dim(groups_corresponding)[2]

        for (i in 1:n_legends)
        {
            group_names <- unique(groups_corresponding[,i])
            group_names_order <- order(group_names) #gets the indices arranged as if they were in alphabetical order
            group_colors <- unique(GroupColorMatrix[,i])
            group_names <- group_names[group_names_order] #arranges them in alphabetical order
            group_colors <- group_colors[group_names_order] #ensures the color order corresponds with the group name order
            x=rep(1,length(group_names))
            pdf(paste(save_directory,'group_color_legend',as.character(i),'.pdf',sep=''),width=5,height=5)
            barplot(x, col=group_colors,border=NA,xlab=NA,ylab=NA,main=ColGroupsScheme[i],ylim=c(0,4),font.main=1)
            legend(x=1,y=4,legend=group_names,col=group_colors,pch = c(15, 15, 15, 15),xpd=NA,bty='n', inset=c(1,1),pt.cex=2)
            dev.off()
         }
     }
}
