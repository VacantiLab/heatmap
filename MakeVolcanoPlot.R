MakeVolcanoPlot <- function(data_location,group_divisions,ColGroupsScheme,genes_to_label=NULL)
{
#data_file: the file (as a string) where the data to be plotted is
  #the items grouped (e.g. patients) are across the columns and the items plotted (e.g. gene expression) are down the rows
  #the rows and columns of the file have headings
#volcano_directory: the directory where the output volcano plot is to be saved
#group_designations_file: the file where the group designations of the columns are stored
  #the file is tab delimited and has the column names across the top
  #each row is a group naming scheme, for example breast cancer patients can be grouped based on the mRNA or protein profiles
    #thus mRNA and protein groupings could be two rows or group naming schemes in this instance
#group_designations_scheme: specifies the row used in the group_designations_file
  #may be 'PAM50' or 'metabolite' for breast cancer patient groupings
#group_divisions: this is a list of 2 lists defined using list(c(),c())
  #each member list contains the names of the groups as they are separated and grouped for analysis in the volcano plot
  #for example, if you are comparing basal tumors to her2, normal like, luminal a, and luminal b, group_divisions <- list(c('basal'),c('her2','normal like','luminal a','luminal b'))
#Note: currently this function only plots -log10(p-value) vs. log2(ratio), where ratio is the median of the fist specified group division over the median of the second specified group division

#Include pertinent libraries
library(ggplot2) #from ggplot2 package, allows the scatter plot to be made

#create the directory 'output' one level above the current directory for storing the heatmap
VolcanoDirectory <- StoreHeatmap()

#Input data
if (!(class(data)=='data.frame')){data <- paste(data_location,'quantities.txt',sep='')}

#Import the data and only keep selected rows if specified
DATA <- OpenDataFile(data,select_rows)

#Can specify a gene to label in red here:
genes_to_label_is_file <- as.logical(grep(genes_to_label,'.txt'))

if (genes_to_label_is_file)
{
  genes_to_label_df <- read.table(file=genes_to_label,head=FALSE,check.names=FALSE,sep='\t')
  genes_to_label_row <- genes_to_label_df[1,]
  n_genes_to_label <- length(genes_to_label_row)
  genes_to_label <- NULL
  for (i in 1:n_genes_to_label) {genes_to_label <- c(genes_to_label,as.character(genes_to_label_row[1,i]))}
}

#store the patient and gene names
patient_names <- colnames(DATA)
gene_names <- rownames(DATA) #store the gene names

#Match the group names to the samples by referencing the group_key
RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme,group_designations_file,group_color_designations_file,select_groups,replicate_scheme)
DATA <- RetrieveGroups_return[[1]]
groups_corresponding <- t(RetrieveGroups_return[[2]])
GroupColorMatrix <- RetrieveGroups_return[[3]]
COLOR_KEY <- RetrieveGroups_return[[4]]

StatTransformByGroup_return <- StatTransformByGroup(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme=ColGroupsScheme)
DATA <- StatTransformByGroup_return[[1]]
groups_corresponding <- StatTransformByGroup_return[[2]]
GroupColorMatrix <- StatTransformByGroup_return[[3]]


#enter the group designations into a copy of the original DATA data frame
DATA_annotated <- DATA
#for each patient, find its group designation and place it in a new row of the DATA_annotated data frame corresponding to the patient's column
for (i in patient_names)
{
  DATA_annotated['group',i] <- GROUP_DESIGNATIONS[group_designations_scheme,i]
}

group_division_1_members <- DATA_annotated['group',] %in% group_divisions[[1]] #find all of the members of group 1 by cross referencing the 'group' row with the first list in the input group_divisions
group_division_2_members <- DATA_annotated['group',] %in% group_divisions[[2]] #find all of the members of group 2 by cross referencing the 'group' row with the second list in the input group_divisions

#for each gene, find the median expression across each group division
for (i in gene_names)
{
  group_division_1_expression <- as.numeric(DATA[i,group_division_1_members]) #get a vector of the expression values across the 1st group
  group_division_2_expression <- as.numeric(DATA[i,group_division_2_members]) #get a vector of the expression values across the 2nd group
  DATA_annotated[i,'ratio'] <- median(group_division_1_expression)/median(group_division_2_expression) #get a ratio of the mean expression
  DATA_annotated[i,'p_value'] <- t.test(group_division_1_expression,group_division_2_expression,alternative='two.sided',mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)$p.value #find the p-value of a t-test, mu is the difference between the two sets according to the null hypothesis
}

#DATA_t_annotated <<- DATA_t_annotated #ggplot looks for variables in global environment, so EXPRESSION must be global

#create columns for the processed variables that are plotted
DATA_annotated$log2ratio <- log2(DATA_annotated$ratio)
DATA_annotated$nlog10p <- -log10(DATA_annotated$p)

#create a new data frame where the rows without numerical values (or gene expressions) are removed
#this removes the row containing the group names
DATA_annotated_less_groups <- DATA_annotated[gene_names,]

#set the plot parameters
XData <- 'log2ratio' #the name of the column of the data frame to be plotted
YData <- 'nlog10p' #the name of the row of the data frame to be plotted
XAxisLimits <- 1.1*c(min(DATA_annotated_less_groups$log2ratio),max(DATA_annotated_less_groups$log2ratio))
YAxisLimits <- 1.1*c(0,max(DATA_annotated_less_groups$nlog10p))
AxisLabels <- c('Log2 Ratio','-Log10(p-value)')
PlotDimensions <- c(4,4)

library(ggplot2)
ScatterPlot <-  ggplot(DATA_annotated_less_groups,aes_string(x=XData,y=YData)) +
  #geom_point(color='black') +

  #################HERE!!!!!!
  #trying to use the subset which requires the plyr package, need to update other packages and restart r to use it - this is to ensure a selected point is visible on top of the others
  geom_point(color='black') + #specifying one color may remove the legend, May be necessary to make this a string when considering groups
  annotate('point',DATA_annotated_less_groups[genes_to_label,'log2ratio'],DATA_annotated_less_groups[genes_to_label,'nlog10p'],col='red') + #points to highlight in red on the volcano plot, x-points first, then y-points
  annotate('text',DATA_annotated_less_groups[genes_to_label,'log2ratio'],DATA_annotated_less_groups[genes_to_label,'nlog10p'],col='#0088B2',label=rownames(DATA_annotated_less_groups[genes_to_label,]),size=0.75) + #points to highlight in red on the volcano plot, x-points first, then y-points
  #geom_errorbar(aes_string(ymin=y_error_low_name,ymax=y_error_high_name)) + #adds error bars to y-variable
  #geom_errorbarh(aes_string(xmin=x_error_low_name,xmax=x_error_high_name)) + #adds error bars to x-variable
  #geom_text(label = rownames(DATA_annotated_less_groups),size=0.75,color='red') + #if you want to label the points
  coord_cartesian(xlim=XAxisLimits,ylim=YAxisLimits,expand = 0) + #expand=0 removes extra space before and after first and last tick marks
  #If you need to specify the axis tick marks, uncomment below and comment out coord_cartesian line above
  #scale_x_continuous(limits=XAxisLimits) +
  #scale_y_continuous(limits=YAxisLimits,breaks=c(-4,-3,-2,-1,0,1,2,3,4)) +
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
  #geom_smooth(method = "lm", se = FALSE, col='red') + #the fitline?
  #annotate("text",x=CorrelationLabelPosition[1],y=CorrelationLabelPosition[2],label=paste('PC = ',pearson_coefficient,' [',pearson_95_confidence[1],', ',pearson_95_confidence[2],']','\n','SC = ',spearman_coefficient,' [',spearman_p_value,']',sep=''),col='red',size=2.5,hjust=0,vjust=1)
  #annotate("text",x=CorrelationLabelPosition[1],y=CorrelationLabelPosition[2],label=paste('PC = ',pearson_coefficient,' [',pearson_95_confidence[1],', ',pearson_95_confidence[2],']',sep=''),col='red',size=3,hjust=0,vjust=1)

print('work!') #There needs to be something in between ggplot call and ggsave call inside a function. This is some sort of bug in ggplot/ggsave

ggsave(paste(volcano_directory,'volcano.pdf',sep=''), width = PlotDimensions[1], height = PlotDimensions[2], dpi = 120)
#rm(DATA_t_annotated) #remove specified global variable so global environment not affected by function. This causes a warning if ColorIndicator is not specified (that's okay)

return(DATA_annotated_less_groups)
}

#############################################################################
#Supporting Functions

#Function to return TRUE if the row does not have any NAs
NoNA <- function(vector)
{
  NoNA <- !is.na(sum(vector))
}
