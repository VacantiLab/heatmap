MakeBoxPlot <- function(data_location,ColGroupsScheme=FALSE,transformation,select_rows='all',select_groups)
# data_file: A text file where the rows are patients and the columns are genes, includes the directory of where the file is located
# box_plot_directory: The location of where the output box plot will be saved, has a "/" at the end
# group_designations_file: A tab delimited text file where each entry corresponds to the group each column in the data_file belongs to
# SelectGenes: A list of the genes, as strings, that are plotted along the x-axis. The default value of 'all' will plot every gene in data_file
# TakeLog2: If TRUE all expression values are log2 transformed, if false values are not log2 transformed
{
    #Include pertinent libraries
    library(ggplot2) #from ggplot2 package
    library(tidyr) #from tidyr package

    #create the directory 'output' one level above the current directory for storing the heatmap
    BoxDirectory <- StoreHeatmap()

    #Input data
    if (is.character(data)){data <- paste(data_location,'quantities.txt',sep='')}
    group_designations_file <- paste(data_location,'group_key.txt',sep='')
    group_color_designations_file <- paste(data_location,'group_color_key.txt',sep='')

    #Import the data and only keep selected rows if specified
    DATA <- OpenDataFile(data,select_rows)

    #Retrieve the corresponding column groupings and keep only those specified
    #Also take the medians if there is a replicate scheme provided (this may not work if the replicate scheme is the only grouping scheme used)
    RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme,group_designations_file,group_color_designations_file,select_groups,replicate_scheme)
    DATA <- RetrieveGroups_return[[1]]
    groups_corresponding <- RetrieveGroups_return[[2]]
    GroupColorMatrix <- RetrieveGroups_return[[3]]

    #Transform the data as specified
    DATA <- transform_data(DATA,transformation)

    #transpose the DATA for making the boxplot if you want
    BoxPlotDataFrame <- as.data.frame(t(BoxPlotDataFrame))

  group_designations <- read.table(file=group_designations_file,head=TRUE,sep='\t',check.names=FALSE) #check.names=FALSE prevents an 'X' from being added to the numeric column names and prevents special characters from being changed
  rownames(group_designations) <- group_designations[,1]
  group_designations[,1] <- NULL

  for (i in 1:length(patients))
  {
    BoxPlotDataFrame[patients[i],'group'] <- as.character(group_designations[GroupingScheme,patients[i]])
  }

  BoxPlotDataFrameLong <- gather_(BoxPlotDataFrame,'condition','measurement',SelectGenes)
  #gather() turns a wide data frame into a long data frame
  #the first argument is the wide data frame
  #the second arguemnt is the name of the column in the long data frame that is composed of column names of the wide data frame
  #the third argument is the name of the column in the long data frame that contains the corresponding values from the wide data frame
  #the last arguemtns are the names of the columns from the wide data frame that are to be included in the long data frame
  #the underscore after the gather function name allows it to take strings as inputs for the column names

  YLabel <- 'Normalized Protein Expression (Log2 Transformed)'
  XLabel <- ''
  TextSize = 26

  FillColors <- c('#E31A1C','#FB9A99','#1F78B4','#A6CEE3','#33A02C')

  #########################MAKE THIS AN INPUT!!!!!!!!!!!!!########################
  #set what is grouped and what is along the x-axis (these can be switched)
  x_var <- 'condition'
  y_var <- 'measurement'
  color_var <- 'group'
  #sets the order of the condition to that specified - this will change depending what the 'condition' is and what the 'group' is
  BoxPlotDataFrameLong$condition <- factor(BoxPlotDataFrameLong$condition,SelectGenes)

  b <- ggplot(BoxPlotDataFrameLong,aes_string(x=x_var, y=y_var)) +
    #stat_boxplot(geom ='errorbar') +
    geom_boxplot(aes_string(fill=color_var),outlier.colour='black',outlier.size=0.5,width=0.4,position=position_dodge(width=0.4),outlier.shape=NA) +
    #geom_hline(yintercept=0, linetype="dashed", color = "red") + #put a line through y=0
    #scale_x_discrete(expand=c(-0.9,0)) +
    #geom_boxplot(aes_string(fill=x_var),outlier.colour='black',outlier.size=0,width=0.7,position=position_dodge(width=0.6)) +
    #scale_y_continuous(breaks=c(0,1,2,3), limits=c(0,3), expand=c(0,0)) +
    #coord_cartesian(ylim=c(0,2)) +
    #stat_summary(position=position_dodge(width=0.5)) +
    #coord_cartesian(ylim=YRangePlot) +
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
    theme(legend.text = element_text(colour="black", size=TextSize)) +
    theme(legend.position=c(0.1,0.85)) +
    scale_fill_manual(values=FillColors) +
    #theme(legend.key.size = unit(0.2, "cm")) +
    labs(x = XLabel) +
    labs(y = YLabel)
  #aes_string() allows the factors to be specified by strings and ensures they are evaluated within the correct environment (aes() causes all sorts of trouble)

  ggsave(paste(box_plot_directory,'boxplot.pdf',sep=''), width = 15, height = 8, dpi = 300, limitsize=FALSE)
  return(b)
}
