assemble_box_plot <- function(DATA_long,FillColors,BoxDirectory,y_bounds)
{
    #set what is grouped and what is along the x-axis (these can be switched, but then may not be compatible with the rest of the MakeBoxPlot function)
    x_var <- 'gene'
    y_var <- 'value'
    color_var <- 'group'

    YLabel <- 'Abundance'
    XLabel <- ''
    TextSize = 8

    pdf_width <- 5
    pdf_height <- 4
    bar_width <- 0.50
    inter_group_spacing <- NULL
    legend_position <- c(0.25,0.9)


    b <- ggplot(DATA_long,aes_string(x=x_var, y=y_var)) +
         geom_boxplot(aes_string(fill=color_var),outlier.colour='black',outlier.size=0.5,width=bar_width,position=position_dodge(width=inter_group_spacing),outlier.shape=NA,lwd=0.2) +
         theme(axis.text.y=element_text(color='black',size=TextSize)) +
         theme(axis.ticks.y=element_line(colour='black',size=0.5)) +
         theme(axis.ticks.x=element_line(colour='black',size=0.5)) +
         theme(axis.text.x=element_text(color='black',size=TextSize,angle=90,hjust=1,vjust=0.5)) +
         theme(axis.title.y=element_text(color='black',vjust=1,size=TextSize)) +
         theme(axis.title.x=element_text(color='black',vjust=0,size=TextSize)) +
         theme(panel.background=element_rect(fill=NA)) +
         theme(axis.line.x=element_line(colour='black',size=0.5,linetype='solid')) +
         theme(axis.line.y=element_line(colour='black',size=0.5,linetype='solid')) +
         theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + #removes gridlines
         theme(legend.title=element_blank()) +
         theme(legend.key=element_rect(fill=NA)) + #No background color in legend
         theme(legend.text = element_text(colour="black", size=(TextSize-2))) +
         theme(legend.position=legend_position) +
         scale_fill_manual(values=FillColors) +
         #theme(legend.key.size = unit(0.2, "cm")) +
         labs(x = XLabel) +
         labs(y = YLabel) +
         coord_cartesian(ylim=y_bounds) #this must be placed inside coord_cartesian() so points outside of the limits are not discarded in calculating medians and IQRs
         #aes_string() allows the factors to be specified by strings and ensures they are evaluated within the correct environment (aes() causes all sorts of trouble)

    ggsave(paste(BoxDirectory,'boxplot.pdf',sep=''), width = pdf_width, height = pdf_height, dpi = 300, limitsize=FALSE)
    return(b)
}
