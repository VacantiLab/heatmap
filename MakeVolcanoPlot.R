MakeVolcanoPlot <- function(data_location,select_groups,ColGroupsScheme,replicate_scheme,genes_to_label=NULL,select_rows=NULL,transformation='none')
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
#select_groups: this is a list of 2 lists defined using list(c(),c()) or it is an array with two group names
  #each member list contains the names of the groups as they are separated and grouped for analysis in the volcano plot
  #for example, if you are comparing basal tumors to her2, normal like, luminal a, and luminal b, select_groups <- list(c('basal'),c('her2','normal like','luminal a','luminal b'))
  #the second entry is the numerator, the first the denominator
#Note: currently this function only plots -log10(p-value) vs. log2(ratio), where ratio is the median of the fist specified group division over the median of the second specified group division
{

    #Stop the program if the replicate scheme is in the ColGroupsScheme
    #Stop the program if more than one ColGroupsScheme or replicate_scheme is specified
    CheckStop(1,parameters=list(ColGroupsScheme,replicate_scheme))
    
    #Include pertinent libraries
    library(ggplot2) #from ggplot2 package, allows the scatter plot to be made

    #create the directory 'output' one level above the current directory for storing the heatmap
    volcano_directory <- StoreHeatmap()

    #Input data
    if (!(class(data)=='data.frame')){data <- paste(data_location,'quantities.txt',sep='')}
    group_designations_file <- paste(data_location,'group_key.txt',sep='')
    group_color_designations_file <- paste(data_location,'group_color_key.txt',sep='')

    #Import the data and only keep selected rows if specified
    #Retrieve the corresponding column groupings and keep only those specified
    #Also take the medians specified by the given replicate scheme
    #There must be a ColGroupsScheme specified, the volcano plot compares two groups
    DATA <- OpenDataFile(data,select_rows)

    #If the groups are packaged into lists, unpack them into a single array
        #This function does not affect the input if it is an array
    UnpackGroups_return <- UnpackGroups(select_groups)
    select_groups <- UnpackGroups_return[[1]]
    group_divisions <- UnpackGroups_return[[2]]

    #Can specify a gene to label in red here:
    genes_to_label_is_file <- FALSE
    if (is.character(genes_to_label))
    {
        genes_to_label_is_file <- length(grep('.txt',genes_to_label)) > 0
    }

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
    #This will also combine replicates if specified to do so, in volcano plots select groups indicates replicates!
    RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme,group_designations_file,group_color_designations_file,select_groups,replicate_scheme)
    DATA <- RetrieveGroups_return[[1]]
    groups_corresponding <- RetrieveGroups_return[[2]]
    GroupColorMatrix <- RetrieveGroups_return[[3]]
    COLOR_KEY <- RetrieveGroups_return[[4]]
    group_values <- RetrieveGroups_return[[5]]

    #If you are concatonating groups, name the new groups and replace all of the groups they map to those with names
    #Also get corresponding colors for those new groups by taking the color that maps to the first sub-group of each concatonated group
    #This function does not affect the input if group_divisions is NULL (i.e. select_groups was not passed as a list to the original function)
    ConcatonateGroups_return <- ConcatonateGroups(group_divisions,groups_corresponding,GroupColorMatrix,COLOR_KEY)
    groups_corresponding <- ConcatonateGroups_return[[1]]
    #GroupColorMatrix <- ConcatonateGroups_return[[2]]

    #Transform the data as specified
    DATA <- transform_data(DATA,transformation)
    gene_name <- rownames(DATA)
    n_gene <- length(gene_name)

    #for each gene, find the median expression across each group division
    for (i in gene_names)
    {
        #DATA[i,'ratio'] <- DATA[i,2]/DATA[i,1] #get a ratio of the mean expression
        group_division_1_expression <- group_values[[1]][i,]
        group_division_2_expression <- group_values[[2]][i,]
        DATA[i,'p_value'] <- t.test(group_division_1_expression,group_division_2_expression,alternative='two.sided',mu=0,paired=FALSE,var.equal=TRUE,conf.level=0.95)$p.value #find the p-value of a t-test, mu is the difference between the two sets according to the null hypothesis
    }

    groups = unique(groups_corresponding)

    #create columns for the processed variables that are plotted
    group1 <- groups[1]
    group2 <- groups[2]
    DATA[,'ratio'] <- DATA[,group2]/DATA[,group1]
    DATA[,'log2ratio'] <- log2(DATA$ratio)
    DATA[,'nlog10p'] <- -log10(DATA$p)

    #set the plot parameters
    XData <- 'log2ratio' #the name of the column of the data frame to be plotted
    YData <- 'nlog10p' #the name of the row of the data frame to be plotted
    XAxisLimits <- 1.1*c(min(DATA$log2ratio),max(DATA$log2ratio))
    YAxisLimits <- 1.1*c(0,max(DATA$nlog10p))
    AxisLabels <- c('Log2 Ratio','-Log10(p-value)')
    PlotDimensions <- c(4,4)

    ScatterPlot <-  ggplot(DATA,aes_string(x=XData,y=YData)) +
        geom_point(color='black') + #specifying one color may remove the legend, May be necessary to make this a string when considering groups
        annotate('point',DATA[genes_to_label,'log2ratio'],DATA[genes_to_label,'nlog10p'],col='red') + #points to highlight in red on the volcano plot, x-points first, then y-points
        annotate('text',DATA[genes_to_label,'log2ratio'],DATA[genes_to_label,'nlog10p'],col='#0088B2',label=rownames(DATA[genes_to_label,]),size=0.75) + #points to highlight in red on the volcano plot, x-points first, then y-points
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

    ggsave(paste(volcano_directory,'volcano.pdf',sep=''), width = PlotDimensions[1], height = PlotDimensions[2], dpi = 120)

    return(DATA)
}

#############################################################################
#Supporting Functions

#Function to return TRUE if the row does not have any NAs
NoNA <- function(vector)
{
    NoNA <- !is.na(sum(vector))
}
