MakeHeatMap <- function(data_location,ColGroupsScheme=FALSE,transformation=NULL,break_seq=NULL,replicate_scheme=NULL,DistanceMethod='pearson',ClusterMethod='ward.D2',data=NULL,select_rows=NULL,select_groups=NULL,label_rows=TRUE,label_cols=TRUE,rev_c_dend=FALSE)
# data_location: a pathway to where the text file containing the data is stored, must have '/' at the end
#    The data file must be named quantities.txt with the genes down the rows and sample names across the columns
#    There must also be a group_key.txt file with the sample names down the rows and the grouping schemes across the columns
#        The group names within each grouping scheme must not match any of those of another grouping scheme
#        Each column must have a header
#    There must also be a group_color_key.txt file with the group names down the rows and the colors beside them
#        All of the group names are listed down the rows with no indication of their scheme membership (thus the names must be unique)
#        Each column must have a header
#    These files must all be tab delimited
# ColGroupsScheme: the name of the grouping scheme used, indicates which row to take from group_key.txt
#    There can only be one for the MakeBoxPlot() function
# transformation: This specifies how the rows should be transformed.
#    Options include: 'log2', 'median_center_iqr_norm', and 'median_norm_log2_transform'
# data: Is a data frame containing the data to be plotted if the data is passed as a data frame instead of through the quantities.txt file in data_location
#    This option could use more testing
# select_rows: These are the genes that you want to plot
#    #If it is NULL, all genes are plotted
# select_groups: This can be an array of the group names (from the provided ColGroupsScheme) that are to be plotted or a list of arrays of group names or NULL
#    If it is an array of group names, those groups are the only ones plotted
#    If it is a list of arrays of group names, groups in the same array are combined into a single group
#    The ratio plotted is the first group median over the second group median
# replicate_scheme: This specifies the grouping scheme that is used to specify groups of replicates
#    This must NOT be a member of ColGroupsScheme, though it must be a grouping scheme defined in group_key.txt
#        As such each member of this grouping scheme must also have colors specified in group_color_key.txt
#    If this is specified, all members of a single group are treated as a single sample and the median values are used
# genes_to_label is either an array of strings corresponding to the genes to be highlighed and labeled in the volcano plot
#    or it is a string that is the name of a text file (ending in .txt) contained within the data_location/gene_lists file that is a column list of the genes to label

{
    #Stop the program if the replicate scheme is in the ColGroupsScheme
    #Stop the program if more than one ColGroupsScheme or replicate_scheme is specified
    CheckStop(1,parameters=list(ColGroupsScheme,replicate_scheme))

    #Include pertinent libraries
    library(ggplot2) #from ggplot2 package, allows the boxplot to be made

    #create the directory 'output' one level above the current directory for storing the heatmap
    output_directory <- StoreHeatmap()

    #Input data
    if (!(class(data)=='data.frame')){data <- paste(data_location,'quantities.txt',sep='')}
    group_designations_file <- paste(data_location,'group_key.txt',sep='')
    group_color_designations_file <- paste(data_location,'group_color_key.txt',sep='')

    #Import the data and only keep selected rows if specified
    DATA <- OpenDataFile(data,select_rows)

    #Retrieve the corresponding column groupings and keep only those specified
    #Also take the medians if there is a replicate scheme provided (this will not work if the replicate scheme is the only grouping scheme used)
    #There must be a ColGroupsScheme specified, as of now it cannot be FALSE

    #If the groups are packaged into lists, unpack them into a single array
        #This function does not affect the input if it is an array
    UnpackGroups_return <- UnpackGroups(select_groups)
    select_groups <- UnpackGroups_return[[1]] #this is now an array of the groups if it were input as a list
    group_divisions <- UnpackGroups_return[[2]] #this is the original input if it was a list, NULL if the input was not a list

    #Match the group names to the samples by referencing the group_key
    ColGroupsScheme_concat <- c(ColGroupsScheme,replicate_scheme)
    RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme_concat,group_designations_file,group_color_designations_file,select_groups)
    groups_corresponding <- RetrieveGroups_return[[1]]
    GroupColorMatrix <- RetrieveGroups_return[[2]]
    COLOR_KEY <- RetrieveGroups_return[[3]]
    CheckStop(5,parameters=list(COLOR_KEY)) #makes sure each group name only has one color assignment

    #See if all of the specified input groups are actually specified in group_key.txt file
    CheckStop(2,parameters=list(select_groups,groups_corresponding))

    #Select the groups that are considered for this box plot
    SelectGroups_return <- SelectGroups(select_groups,DATA,ColGroupsScheme_concat,groups_corresponding,GroupColorMatrix,inclusion_grouping_scheme=ColGroupsScheme)
    #inclusion_grouping_scheme will need to be specified when more than one grouping scheme can be used such as in a heatmap
    DATA <- SelectGroups_return[[1]]
    groups_corresponding <- SelectGroups_return[[2]]
    GroupColorMatrix <- SelectGroups_return[[3]]

    #If you are concatonating groups, name the new groups and replace all of the groups they map to those with names
    #Also get corresponding colors for those new groups by taking the color that maps to the first sub-group of each concatonated group
    #This function does not affect the input if group_divisions is NULL (i.e. select_groups was not passed as a list to the original function)
    ConcatonateGroups_return <- ConcatonateGroups(group_divisions,groups_corresponding,GroupColorMatrix,COLOR_KEY,concat_group_scheme=ColGroupsScheme)
    #concat_group_scheme has to be ColGroupsScheme because a volcano plot only has one ColGroupsScheme
    groups_corresponding <- ConcatonateGroups_return[[1]]
    GroupColorMatrix <- ConcatonateGroups_return[[2]]
    groups_concatonated <- ConcatonateGroups_return[[3]]
    colors_concatonated <- ConcatonateGroups_return[[4]]
    group_concationation <- is.list(group_divisions)

    #return median of replicates if specified to do so
    MedianGroup_return <- MedianGroup(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme,ColGroupsScheme)
    DATA <- MedianGroup_return[[1]]
    groups_corresponding <- MedianGroup_return[[2]]
    GroupColorMatrix <- MedianGroup_return[[3]]

    #Transform the data as specified
    DATA <- transform_data(DATA,transformation)
    gene_name <- rownames(DATA)
    n_gene <- length(gene_name)

    #Retrieve the clusters for creating the dendrograms
    ClusterData_return <- ClusterData(DATA,DistanceMethod,ClusterMethod,rev_c_dend)
    colv <- ClusterData_return[[1]]
    rowv <- ClusterData_return[[2]]
    C_col <- ClusterData_return[[3]]
    C_row <- ClusterData_return[[4]]

    #make numeric matrix out of data frame
    DifExpMatx=matrix(as.numeric(unlist(DATA)),nrow=nrow(DATA))
    rownames(DifExpMatx)<-rownames(DATA) #name the rows of the numeric matrix so they are transfered to the dendrogram in the heat map
    colnames(DifExpMatx)<-colnames(DATA) #name the columns of the numeric matrix so they are transfered to the dendrogram in the heat map

    #Get the row names and column names for the heatmap
    label_rows <- FindRowLabels(label_rows,DifExpMatx)
    label_cols <- FindColLabels(label_cols,DifExpMatx)

    #Draw and save the heatmap
    assemble_heatmap_return <- assemble_heatmap(GroupColorMatrix,DifExpMatx,colv,rowv,heat_map_colors,break_seq,label_rows,label_cols,output_directory,DistanceMethod,ClusterMethod)
    heat_map_colors <- assemble_heatmap_return[[1]]
    dev.off() #turn off printing to the specified pdf

    #Draw and save the legends
    MakeColorKey(break_seq,heat_map_colors,output_directory)
    MakeGroupLegend(groups_corresponding,GroupColorMatrix,ColGroupsScheme,output_directory)

    #Return what could be used
    return(list(C_col,C_row,groups_corresponding,DATA))
}
