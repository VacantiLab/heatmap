assemble_volcano_plot <- function(rp_df,output_directory)
# rp_df is ratio p-value data frame
{
#set the plot parameters
XData <- 'log2_ratio' #the name of the column of the data frame to be plotted
YData <- 'nlog10_p' #the name of the row of the data frame to be plotted
XAxisLimits <- 1.1*c(min(rp_df$log2_ratio),max(rp_df$log2_ratio))
YAxisLimits <- 1.1*c(0,max(rp_df$nlog10_p))
AxisLabels <- c('Log2 Ratio','-Log10(p-value)')
PlotDimensions <- c(4,4)

ScatterPlot <-  ggplot(rp_df,aes_string(x=XData,y=YData)) +
    geom_point(color='black') + #specifying one color may remove the legend, May be necessary to make this a string when considering groups
    #annotate('point',rp_df[genes_to_label,XData],rp_df[genes_to_label,YData],col='red') + #points to highlight in red on the volcano plot, x-points first, then y-points
    #annotate('text',rp_df[genes_to_label,XData],rp_df[genes_to_label,YData],col='#0088B2',label=rownames(DATA[genes_to_label,]),size=0.75) + #points to highlight in red on the volcano plot, x-points first, then y-points
    coord_cartesian(xlim=XAxisLimits,ylim=YAxisLimits,expand = 0) + #expand=0 removes extra space before and after first and last tick marks
    theme(axis.text.y=element_text(color='black',size=10)) +
    theme(axis.ticks.y=element_line(colour='black',size=1)) +
    theme(axis.ticks.x=element_line(colour='black',size=1)) +
    theme(axis.text.x=element_text(color='black',size=10)) +
    theme(axis.title.y=element_text(color='black',vjust=1,size=10)) +
    theme(axis.title.x=element_text(color='black',vjust=0,size=10)) +
    theme(panel.background=element_rect(fill=NA)) +
    theme(axis.line.x=element_line(colour='black',size=1,linetype='solid')) +
    theme(axis.line.y=element_line(colour='black',size=1,linetype='solid')) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +#removes gridlines
    theme(legend.title=element_blank()) +
    theme(legend.key=element_rect(fill=NA)) +
    theme(legend.text = element_text(colour="black", size = 6)) +
    #theme(legend.position=LegendPosition) +
    theme(legend.key.size = unit(0.3, "cm")) + #Space between legend symbols and text, maybe?
    theme(legend.background = element_rect(fill="transparent",linetype = 0,colour = 'transparent')) +
    labs(x=AxisLabels[1]) +
    labs(y=AxisLabels[2])

print('work!') #There needs to be something in between ggplot call and ggsave call inside a function. This is some sort of bug in ggplot/ggsave

ggsave(paste(output_directory,'volcano.pdf',sep=''), width = PlotDimensions[1], height = PlotDimensions[2], dpi = 120)

return(ScatterPlot)
}
