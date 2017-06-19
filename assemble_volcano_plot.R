assemble_volcano_plot <- function(rp_df,output_directory,genes_to_label,XData,YData,filename)
# rp_df is ratio p-value data frame
{
#set the plot parameters
XAxisLimits <- 1.1*c(min(rp_df[,XData]),max(rp_df[,XData]))
if (min(rp_df[,YData])>0){ymin=0}
if (min(rp_df[,YData])<0){ymin=1.1*min(rp_df[,YData])}
YAxisLimits <- 1.1*c(ymin,max(rp_df[,YData]))
AxisLabels <- c('Log2 Ratio','-Log10(p-value)')
PlotDimensions <- c(4,4)

#linear_fit <- geom_smooth(method = "lin_model", se = FALSE, col='red',size=0.1) +
if (filename == 'volcano.pdf'){linear_fit <- NULL}
if (filename == 'volcano_regression.pdf'){linear_fit <- geom_smooth(method = "lm", se = FALSE, col='red',size=0.1)} #lm is a function used to fit a line

ScatterPlot <-  ggplot(rp_df,aes_string(x=XData,y=YData)) +
    geom_point(color='black') + #specifying one color may remove the legend, May be necessary to make this a string when considering groups
    annotate('point',rp_df[genes_to_label,XData],rp_df[genes_to_label,YData],col='red') + #points to highlight in red on the volcano plot, x-points first, then y-points
    annotate('text',rp_df[genes_to_label,XData],rp_df[genes_to_label,YData],col='#0088B2',label=rownames(rp_df[genes_to_label,]),size=0.75) + #points to highlight in red on the volcano plot, x-points first, then y-points
    linear_fit +
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

ggsave(paste(output_directory,filename,sep=''), width = PlotDimensions[1], height = PlotDimensions[2], dpi = 120)

return(ScatterPlot)
}