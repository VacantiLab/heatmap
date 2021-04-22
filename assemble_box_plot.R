assemble_box_plot <- function(DATA_long,
                              FillColors,
                              output_directory,
                              ybounds,
                              qc_plot,
                              box_plot_type,
                              plot_width,
                              plot_height,
                              bar_width,
                              legend_position,
                              text_angle,
                              transformation,
                              ytick,
                              ErrorBarSize=0.75,
                              PointSize=3,
                              ErrorFile,
                              x_var,
                              y_var,
                              color_var)
{
    #plot_type: can be 'boxplot', 'scatter_bar_plot', or 'bar_plot'


    if (qc_plot)
    {
        x_var <- 'sample'
    }

#    x_var <- 'group'
#    y_var <- 'value'
#    color_var <- 'gene'

    no_groups_color <- 'grey'

    inter_group_spacing <- 1.1*bar_width

    YLabel <- ''
    XLabel <- ''
    TextSize = 10

    pdf_width <- unit(plot_width,'in')
    pdf_height <- unit(plot_height,'in')

    #If there are groups, the group color is specified by the color_var
    if(!is.null(FillColors))
    {
        #box_plot_type='boxplot'
        gbp <- geom_boxplot(aes_string(fill=color_var),
                            outlier.colour='black',
                            outlier.size=0.5,
                            width=bar_width,
                            position=position_dodge(width=inter_group_spacing),
                            outlier.shape=NA,
                            lwd=0.2)

        #If you want labeled outliers
        #gbp <- geom_boxplot(aes_string(fill=color_var),outlier.colour='black',outlier.size=0.5,width=bar_width,position=position_dodge(width=inter_group_spacing),outlier.shape=20,lwd=0.2)

        #need this for point-errorbar format
        #If you want a scatter plot with mean and SD
        library(plyr) #this package contains the ddply function which allows for making a data frame with summary statistics

        # If a separate file is specified with error values, use them
        if (is.character(ErrorFile))
        {
            ERROR <- read_txt_to_df(ErrorFile)
            for(row in 1:length(DATA_long[,'value']))
            {
              gene <- DATA_long[row,'gene']
              group <- DATA_long[row,'group']
              DATA_long[row,'sd'] <- ERROR[gene,group]
            }

            DATA_long_summary <- ddply(DATA_long,
                                      c(color_var,x_var),
                                      summarise,
                                      value2=mean(value),
                                      sd=sd)
        }

        # If no separate file is specified with error values, calculate the error values
        if (!is.character(ErrorFile))
        {
            DATA_long_summary <- ddply(DATA_long,
                                      c(color_var,x_var),
                                      summarise,
                                      value2=mean(value),
                                      sd=sd(value))
        }


        #  if you use 'value' to name the column instead of 'value2', the standard deviations will not calculate
        colnames(DATA_long_summary)[colnames(DATA_long_summary)=='value2']='value'
        #  putting the name, 'value', back as the column name

        #box_plot_type='scatter_bar_plot'
        gep <- geom_errorbar(aes(ymin=value-sd,ymax=value+sd),width=0.75*bar_width,position=position_dodge(width=inter_group_spacing))
        gpp <- geom_point(aes(color=group),position=position_dodge(width=inter_group_spacing),size=PointSize)

        #box_plot_type='bar_plot'
        grp <- geom_bar(aes(fill=group),position=position_dodge(width=inter_group_spacing),width=bar_width,stat='identity')

    }

    #If there are no groups, the group color is specified by no_groups_color
    if(is.null(FillColors))
    {
        gbp <- geom_boxplot(outlier.colour='black',outlier.size=0.5,width=bar_width,position=position_dodge(width=inter_group_spacing),outlier.shape=NA,lwd=0.2,fill=no_groups_color)
        #gbp <- geom_boxplot(outlier.colour='black',outlier.size=0.5,width=bar_width,position=position_dodge(width=inter_group_spacing),outlier.shape=20,lwd=0.2,fill=no_groups_color)
        # second line provides outlier shape
    }

    #If you want to add markers for points on the boxplot
    boxplot_points <- FALSE #for now this needs to be changed in the source code to have points on the boxplot
    gtp <- NULL
    if(boxplot_points)
    {
        #If you want to include all points
        gtp <- geom_text(aes_string(label='sample',group=color_var),size=2,position=position_dodge(width=inter_group_spacing))

        #If you want the point to be the mean - the points are colored
        #gtp <- geom_point(aes(color=group),position=position_dodge(width=inter_group_spacing),size=2.5)
    }

    axis_limits <- coord_cartesian(ylim=ybounds) #this must be placed inside coord_cartesian() so points outside of the limits are not discarded in calculating medians and IQRs
    y_scale <- NULL # this is a graph property defined below, not the same as ytick, the tick mark definition
    text_angle_indicator <- NULL

    if(box_plot_type=='boxplot'){gpp<-NULL; gep<-NULL; grp<-NULL; DATA_to_plot <- DATA_long}
    if(box_plot_type=='scatter_bar_plot')
    {
        grp<-NULL
        gbp<-NULL
        DATA_to_plot <- DATA_long_summary
        if(!is.null(FillColors))
        {
            gep <- geom_errorbar(aes(ymin=value-sd,ymax=value+sd,color=group),width=0.75*bar_width,position=position_dodge(width=inter_group_spacing),size=ErrorBarSize)
        }

    }
    if(box_plot_type=='bar_plot')
    {
        gpp<-NULL
        gbp<-NULL
        DATA_to_plot <- DATA_long_summary

        upper_y_bound <- ceiling(ybounds[2]*10)*0.1 #the ybound as it's calculated from the data or specified and then rounded up to an even tenth
        if (is.null(ytick)){ytick = upper_y_bound/2}
        if (!grepl('log',transformation) || is.null(transformation)){ybounds <- c(0,upper_y_bound)} #only set 0 to ymin if not log transformed
        x_bounds <- c(0.5,length(unique(DATA_long$gene))+0.5) #when used with expand=F in coord_cartesian, sets a little space on either side of x-axis variables
        axis_limits <- coord_cartesian(ylim=ybounds, expand=F, xlim=x_bounds) #this must be placed inside coord_cartesian() so points outside of the limits are not discarded in calculating medians and IQRs
        y_scale <- scale_y_continuous(breaks = seq(ybounds[1],ybounds[2], by=ytick)) #scale_y_continuous can be used with coord_cartesian without affecting data calculations
    }

    #axis label placement
    hjust_value <- 0
    if (text_angle==90 | text_angle==45){hjust_value <- 1; vjust_value <- 1}
    if (text_angle==0){hjust_value <- 0.5; vjust_value <- 0.5}
    text_angle_indicator <- theme(axis.text.x = element_text(angle=text_angle,hjust=hjust_value,vjust=vjust_value))
    
    if (color_var == 'group')
    {
        ScaleFillDesignation <- scale_fill_manual(values=FillColors)
        ScaleColorDesignation <- scale_color_manual(values=FillColors)
    }
    
    if (color_var == 'gene')
    {
        ScaleFillDesignation <- NULL
        ScaleColorDesignation <- NULL
    }

    b <- ggplot(DATA_to_plot,aes_string(x=x_var, y=y_var, fill=color_var)) +
         gpp +
         gbp +
         gtp +
         grp +
         gep +
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
         theme(legend.text = element_text(colour="black", size=(TextSize))) +
         theme(legend.position=legend_position) +
         theme(plot.margin = margin(10,10,0,0)) +
         ScaleFillDesignation +
         ScaleColorDesignation +
         #theme(legend.key.size = unit(0.2, "cm")) +
         labs(x = XLabel) +
         labs(y = YLabel) +
         axis_limits + #this must be placed inside coord_cartesian() so points outside of the limits are not discarded in calculating medians and IQRs
         y_scale +
         text_angle_indicator
         #aes_string() allows the factors to be specified by strings and ensures they are evaluated within the correct environment (aes() causes all sorts of trouble)

    ggsave(paste(output_directory,'boxplot.pdf',sep=''), width = pdf_width, height = pdf_height, dpi = 300, limitsize=FALSE)
    return(b)
}
