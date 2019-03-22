assemble_box_plot <- function(DATA_long,FillColors,output_directory,y_bounds,qc_plot)
{
    #set what is grouped and what is along the x-axis (these can be switched, but then may not be compatible with the rest of the MakeBoxPlot function)
    x_var <- 'gene'
    y_var <- 'value'
    color_var <- 'group'

    if (qc_plot)
    {
        x_var <- 'sample'
    }

#    x_var <- 'group'
#    y_var <- 'value'
#    color_var <- 'gene'

    no_groups_color <- 'grey'

    YLabel <- 'Relative Abundance'
    XLabel <- ''
    TextSize = 8

    pdf_width <- unit(3,'in')
    pdf_height <- unit(2,'in')
    bar_width <- 0.30
    inter_group_spacing <- 0.40
    legend_position <- c(0.10,0.92)

    #If there are groups, the group color is specified by the color_var
    if(!is.null(FillColors))
    {
        gbp <- geom_boxplot(aes_string(fill=color_var),outlier.colour='black',outlier.size=0.5,width=bar_width,position=position_dodge(width=inter_group_spacing),outlier.shape=NA,lwd=0.2)

        #If you want labeled outliers
        #gbp <- geom_boxplot(aes_string(fill=color_var),outlier.colour='black',outlier.size=0.5,width=bar_width,position=position_dodge(width=inter_group_spacing),outlier.shape=20,lwd=0.2)

        ##need this for point-errorbar format
        ##If you want a scatter plot with mean and SD
        #library(plyr) #this package contains the ddply function which allows for making a data frame with summary statistics
        #DATA_long_summary <- ddply(DATA_long,c(color_var,x_var),summarise,value2=mean(value),sd=sd(value))
        ##  if you use 'value' to name the column instead of 'value2', the standard deviations will not calculate
        #colnames(DATA_long_summary)[colnames(DATA_long_summary)=='value2']='value'
        ##  putting the name, 'value', back as the column name
        #gbp <- geom_errorbar(aes(ymin=value-sd,ymax=value+sd,color=group),width=bar_width,position=position_dodge(width=inter_group_spacing))

    }

    #If there are no groups, the group color is specified by no_groups_color
    if(is.null(FillColors))
    {
        gbp <- geom_boxplot(outlier.colour='black',outlier.size=0.5,width=bar_width,position=position_dodge(width=inter_group_spacing),outlier.shape=NA,lwd=0.2,fill=no_groups_color)
        #gbp <- geom_boxplot(outlier.colour='black',outlier.size=0.5,width=bar_width,position=position_dodge(width=inter_group_spacing),outlier.shape=20,lwd=0.2,fill=no_groups_color)
        # second line provides outlier shape
    }

    #If you want to add markers for points on the boxplot
    boxplot_points <- TRUE #for now this needs to be changed in the source code to have points on the boxplot
    gtp <- NULL
    if(boxplot_points)
    {
        #If you want to include all points
        #gtp <- geom_text(aes_string(label='sample',group=color_var),size=2,position=position_dodge(width=inter_group_spacing))

        #If you want the point to be the mean - the points are colored
        gtp <- geom_point(aes(color=group),position=position_dodge(width=inter_group_spacing),size=2.5)
    }

    #set y-bounds manually if desired
    manual_ybounds <- FALSE
    if(manual_ybounds){y_bounds <- c(-2,3)}

    #need this for point-errorbar format
    #scatter_box here indicates that you just want points and error bars
    #    for now this is hardcoded
    scatter_box <- FALSE
    if(scatter_box){DATA_to_plot <- DATA_long_summary}
    if(!scatter_box){DATA_to_plot <- DATA_long}

    b <- ggplot(DATA_to_plot,aes_string(x=x_var, y=y_var)) +
         #gvp +
         gbp +
         gtp +
         theme(axis.text.y=element_text(color='black',size=TextSize)) +
         theme(axis.ticks.y=element_line(colour='black',size=0.5)) +
         theme(axis.ticks.x=element_line(colour='black',size=0.5)) +
         #Set the x-axis text
         #theme(axis.text.x=element_text(color='black',size=TextSize,angle=90,hjust=1,vjust=0.5)) +
         theme(axis.text.x=element_text(color='black',size=TextSize,angle=0,vjust=0.5)) +
         theme(axis.title.y=element_text(color='black',vjust=1,size=TextSize)) +
         theme(axis.title.x=element_text(color='black',vjust=0,size=TextSize)) +
         theme(panel.background=element_rect(fill=NA)) +
         theme(axis.line.x=element_line(colour='black',size=0.5,linetype='solid')) +
         theme(axis.line.y=element_line(colour='black',size=0.5,linetype='solid')) +
         theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + #removes gridlines
         theme(legend.title=element_blank()) +
         theme(legend.background=element_rect(fill=NA)) + #No background color in legend
         theme(legend.key.size=unit(0.3,'cm')) +
         theme(legend.text = element_text(colour="black", size=(TextSize-2))) +
         theme(legend.position=legend_position) +
         scale_fill_manual(values=FillColors) +
         scale_color_manual(values=FillColors) +
         #theme(legend.key.size = unit(0.2, "cm")) +
         labs(x = XLabel) +
         labs(y = YLabel) +
         coord_cartesian(ylim=y_bounds) #this must be placed inside coord_cartesian() so points outside of the limits are not discarded in calculating medians and IQRs
         #aes_string() allows the factors to be specified by strings and ensures they are evaluated within the correct environment (aes() causes all sorts of trouble)

    ggsave(paste(output_directory,'boxplot.pdf',sep=''), width = pdf_width, height = pdf_height, dpi = 300, limitsize=FALSE)
    return(b)
}
