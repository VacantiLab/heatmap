assemble_scatter_plot <- function(DATA,FillColors,output_directory,select_rows)
{
    #set the gene names across columns and the sample names down the rownames
    DATA_t <- as.data.frame(t(DATA))

    #set what is grouped and what is along the x-axis (these can be switched, but then may not be compatible with the rest of the MakeBoxPlot function)
    x_var <- select_rows[1]
    y_var <- select_rows[2]
    color_var <- 'group'
    no_groups_color <- 'grey'

    XLabel <- select_rows[1]
    YLabel <- select_rows[2]
    TextSize = 8

    pdf_width <- unit(4,'cm')
    pdf_height <- unit(4,'cm')

    #If there are groups, the group color is specified by the color_var, right now group colors cannot be specified
    if(!is.null(FillColors))
    {
        gsp <- geom_point(color=no_groups_color)
    }

    #If there are no groups, the group color is specified by no_groups_color
    if(is.null(FillColors))
    {
        gsp <- geom_point(color=no_groups_color)
    }

    ScatterPlot <-  ggplot(DATA_t,aes_string(x=x_var,y=y_var)) +
      gsp +
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
      theme(legend.key.size = unit(0.3, "cm")) + #Space between legend symbols and text, maybe?
      theme(legend.background = element_rect(fill="transparent",linetype = 0,colour = 'transparent'))

   ggsave(paste(output_directory,'scatterplot.pdf',sep=''), width = pdf_width, height = pdf_height, dpi = 300, limitsize=FALSE)
}
