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
    #Extract the data required to make a heatmap
    ArrangeData_return <- ArrangeData(ColGroupsScheme,replicate_scheme,transformation,data,data_location,select_rows,select_groups,visualization='heatmap')
    sig_test_list <- ArrangeData_return[[1]]
    output_directory <- ArrangeData_return[[2]]
    group_order <- ArrangeData_return[[3]]
    gene_name <- ArrangeData_return[[4]]
    DATA_long <- ArrangeData_return[[5]]
    FillColors <- ArrangeData_return[[6]]
    DATA <- ArrangeData_return[[7]]
    GroupColorMatrix <- ArrangeData_return[[8]]
    groups_corresponding <- ArrangeData_return[[9]]

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
    assemble_heatmap_return <- assemble_heatmap(GroupColorMatrix,DifExpMatx,colv,rowv,break_seq,label_rows,label_cols,output_directory,DistanceMethod,ClusterMethod)
    heat_map_colors <- assemble_heatmap_return[[1]]
    dev.off() #turn off printing to the specified pdf

    #Draw and save the legends
    MakeColorKey(break_seq,heat_map_colors,output_directory)
    MakeGroupLegend(groups_corresponding,GroupColorMatrix,ColGroupsScheme,output_directory)

    #Return what could be used
    return(list(C_col,C_row,groups_corresponding,DATA))
}
