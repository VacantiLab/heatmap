MakeHeatMap <- function(dl,ColGroupsScheme=NULL,transformation='log2',break_seq=seq(-2,2,0.5),replicate_scheme=NULL,DistanceMethod='pearson',ClusterMethod='ward.D2',data=NULL,select_rows=NULL,select_groups=NULL,inclusion_grouping_scheme=NULL,label_rows=FALSE,label_cols=FALSE,rev_c_dend=FALSE,ddt=NULL,med_norm=FALSE,handle_blanks='remove_row',presentation='normal',visualization='heatmap')
# DL: stands for data location, a pathway to where the text file containing the data is stored, must have '/' at the end
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
#    If it is a path to a .txt file, that file must contain the gene names, one on each line, and they will be used as the select_rows
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
# ddt: data dependent transform
#    Allows for the transformation of the columns to ratios of columns
#    It is a list of two arrays
#        The first array contains the names of the columns of the original DATA set that will be the numerators
#        The second array contains the names of the columns of the original DATA set that will be the denominators
#    It can also be a string that is one of the members of ColGroupsScheme
#        If this is the case, the data within groups designated by the ddt grouping scheme will be normalized to the median of the group
#        ddt would not actually be a group that is presented in the plot but must be a member of the ColGroupsScheme input
# med_norm specifies to median normalize columns, will occur before rows are selected or other transformations performed (Done within OpenDataFile)
# handle_blanks is used in OpenDataFile in the SupportingFunctions.R file. It is a string that specifies what to do with blank values
# presentation: can be set to 'correlation_matrix' to display the correlation of each row vs. each other row as the matrix visualization
#    The original clustering of the rows of the input data will determine the visualization dendrograms (along the rows and colums - they're the same)

{
    #Extract the data required to make a heatmap
    data_location <- dl
    ArrangeData_return <- ArrangeData(ColGroupsScheme,replicate_scheme,transformation,data,data_location,select_rows,select_groups,visualization=visualization,ddt,med_norm,handle_blanks,inclusion_grouping_scheme,ttest=FALSE)
    sig_test_list <- ArrangeData_return[[1]]
    output_directory <- ArrangeData_return[[2]]
    group_order <- ArrangeData_return[[3]]
    gene_name <- ArrangeData_return[[4]]
    DATA_long <- ArrangeData_return[[5]]
    FillColors <- ArrangeData_return[[6]]
    DATA <- ArrangeData_return[[7]]
    GroupColorMatrix <- ArrangeData_return[[8]]
    groups_corresponding <- ArrangeData_return[[9]]
    DATA_original <- ArrangeData_return[[11]] #the data before any transformations or exclusions were made


    #Retrieve the clusters for creating the dendrograms
    ClusterData_return <- ClusterData(DATA,DistanceMethod,ClusterMethod,rev_c_dend)
    colv <- ClusterData_return[[1]]
    rowv <- ClusterData_return[[2]] #dendrogram object
    C_col <- ClusterData_return[[3]]
    C_row <- ClusterData_return[[4]] #hclust object
    Cor_col <- ClusterData_return[[5]] #column correlation matrix
    Cor_row <- ClusterData_return[[6]] #row correlation matrix

    #Create the edge table if specified to do so
    if(visualization=='edge_table')
    {
        CreateEdgeTable(Cor_row,DATA)
    }

    print('making the heatmap')
    #make numeric matrix out of data frame
    DifExpMatx=matrix(as.numeric(unlist(DATA)),nrow=nrow(DATA))
    rownames(DifExpMatx)<-rownames(DATA) #name the rows of the numeric matrix so they are transfered to the dendrogram in the heat map
    colnames(DifExpMatx)<-colnames(DATA) #name the columns of the numeric matrix so they are transfered to the dendrogram in the heat map

    #Get the row names and column names for the heatmap
    label_rows <- FindRowLabels(label_rows,DifExpMatx)
    label_cols <- FindColLabels(label_cols,DifExpMatx)

    #Draw and save the heatmap
    n_clusters = 7
    assemble_heatmap_return <- assemble_heatmap(GroupColorMatrix,DifExpMatx,colv,rowv,break_seq,label_rows,label_cols,output_directory,DistanceMethod,ClusterMethod,C_col,C_row,Cor_col,Cor_row,presentation,n_clusters)
    heat_map_colors <- assemble_heatmap_return[[1]]
    cluster_colors <- assemble_heatmap_return[[2]]
    cutree_genes <- assemble_heatmap_return[[3]]

    #Draw and save the legends
    MakeColorKey(break_seq,heat_map_colors,output_directory)
    MakeGroupLegend(groups_corresponding,GroupColorMatrix,ColGroupsScheme,output_directory)
    MakeClusterLegend(n_clusters,cluster_colors,output_directory)

    #Return what could be used
    return(list(C_col,C_row,groups_corresponding,DATA,cutree_genes,DATA_original))
}
