MakeHeatMap <- function(dl,
                        ColGroupsScheme=NULL,
                        transformation='log2',
                        break_seq=seq(-2,2,0.5),
                        replicate_scheme=NULL,
                        DistanceMethod='pearson',
                        ClusterMethod='ward.D2',
                        data=NULL,
                        select_rows=NULL,
                        select_groups=NULL,
                        inclusion_grouping_scheme=NULL,
                        label_rows=FALSE,
                        label_cols=FALSE,
                        rev_c_dend=FALSE,
                        ddt=NULL,
                        handle_blanks='remove_row',
                        presentation='normal',
                        visualization='heatmap',
                        select_rows_after_transform=NULL,
                        transform_after_column_exclusion=FALSE,
                        graphics_type='.png',
                        n_clusters=1,
                        FilterRowsForMeanSpread=NULL)
# Inputs:
# dl: stands for data location, a pathway to where the text file containing the data is stored, must have '/' at the end
#    The data file must be named quantities.txt with the genes down the rows and sample names across the columns
#    There must also be a group_key.txt file with the sample names down the rows and the grouping schemes across the columns
#        The group names within each grouping scheme must not match any of those of another grouping scheme
#        Each column must have a header
#    There must also be a group_color_key.txt file with the group names down the rows and the colors beside them
#        All of the group names are listed down the rows with no indication of their scheme membership (thus the names must be unique)
#        Each column must have a header
#    These files must all be tab delimited
# ColGroupsScheme: the name of the grouping scheme used indicates which column to take from group_key.txt
#    There can only be one for the MakeBoxPlot() function
# transformation: This specifies how the rows should be transformed.
#    Options include: 'log2', 'median_center_iqr_norm', and 'median_norm_log2_transform', and others.
#    See transform_data.R for a complete list of possible transformations.
# data: Is a data frame containing the data to be plotted if the data is passed as a data frame
#       If this is used then the data is not passed through the quantities.txt file in data_location
# select_rows: These are the rows that you want to select before transforming - the transformation will be done on these rows
#    If it is NULL, all genes included
#    If it is a path to a .txt file, that file must contain the gene names, one on each line, and they will be used as the select_rows
#    The transformations will not be done on the whole data set and this does impact the values
# select_rows_after_transform: these rows are selected to be plotted after tranformations have been completed
#     This does not impact values
# transform_after_column_exclusion: if TRUE, transformations happen after groups (columns) are selected, not on the whole data set
#     If TRUE, select_groups impacts the values
#     If FALSE, select_groups does not impact the values
#     If TRUE, this may impact values if ddt is specified
#     If FALSE, this will not impact values if ddt is specified
# select_groups: This can be an array of the group names that are to be plotted or a list of arrays of group names or NULL
#    If it is an array of group names, those groups are the only ones plotted
#    If it is a list of arrays of group names, groups in the same array are combined into a single group
#    If it is NULL, all groups in the ColGroupsScheme are plotted
#    These can be groups outside of ColGroupsScheme, but the scheme must then be specified as the inclusion_grouping_scheme
#    The spcification of transform_after_column_exclusion determines if this impacts data transformation/normalization
# inclusion_grouping_scheme: This is the grouping scheme that you want to specify to select the columns with
#    It can remain as the default NULL and the groups will be selected based on the first grouping scheme in ColGroupsScheme if select_groups is specified
#    It can be outside of the ColGroupsScheme as well
# ddt: data dependent transformation; this is a grouping scheme that all samples within that group are normalized to its median and then log2 transformed
#    It must be one of the groups specified in the ColGroupsScheme
#    It is not presented as a ColGroupsScheme, it is just used for normalization purposes
#        For example, say you have cell lines control and treated
#            You can specify to normalize within cell lines and then use select_groups and inclusion_grouping_scheme to plot only the treated samples
#                The result would be the treatment response for each cell line
#    transform_after_column_exclusion must be FALSE because the transformation would also occur after DDT which doesn't make sense
#        Sample loading should be accounted for before DDT
#        The transformation should also be linear because DDT log2 transforms resulting ratios
# replicate_scheme: This specifies the grouping scheme that is used to specify groups of replicates
#    This must NOT be a member of ColGroupsScheme, though it must be a grouping scheme defined in group_key.txt
#        As such each member of this grouping scheme must also have colors specified in group_color_key.txt
#    If this is specified, all members of a single group are treated as a single sample and the median values are used
# ddt: data dependent transform
#    Allows for the transformation of the columns to ratios of columns
#    It is a list of two arrays
#        The first array contains the names of the columns of the original DATA set that will be the numerators
#        The second array contains the names of the columns of the original DATA set that will be the denominators
#    It can also be a string that is one of the members of ColGroupsScheme
#        If this is the case, the data within groups designated by the ddt grouping scheme will be normalized to the median of the group
#        ddt would not actually be a group that is presented in the plot but must be a member of the ColGroupsScheme input
# handle_blanks: Used in OpenDataFile in the SupportingFunctions.R file.
#                It is a string that specifies what to do when blank values are encountered
#                Can be 'remove_row' or 'replace_with_rowmin'
# presentation: can be set to 'correlation_matrix' to display the correlation of each row vs. each other row as the matrix visualization
#    The original clustering of the rows of the input data will determine the visualization dendrograms (along the rows and columns - they're the same)
#    Right now there is a problem with col_side_colors and it doesn't work
# label_rows: TRUE or FALSE to label all rows, a vector to label rows in the vector, or a text file to label rows stated on each line
# graphics_type: can be either '.pdf', '.png', or '.jpeg' and dictates the type of file the heatmap will be printed to
# n_clusters: Specifies the number of clusters on the rows that are identified by different colors
# FilterRowsForMeanSpread: If it is a decimal fraction (such as 0.7), rows are filtered out who have more than 70% of their values below or above 0
#     This is meant for log2 transformed data where there may be a small cluster of columns driving the profile
#     Corrects red or blue streaks across the heat map

# Outputs are a list: MHR
# 1. ADL: "Arrange DATA Return"; A list of items returned by the function ArrangeData()
#     01. sig_test_list:
#     02. output_directory:
#     03. group_order: column dendrogram object
#     04. gene_name: row dendrogram object
#     05. DATA_long:
#     06. FillColors:
#     07. DATA:
#     08. GroupColorMatrix:
#     09. groups_corresponding:
#     10. DATA_original:

# 2. AHR: "Assemble Heatmap Return"; A list of items returned by the function assemble_heatmap()
#     01. heat_map_colors:
#     02. cutree_genes:

# 3. AHI: "Assmble Heatmap Input"; A list of everything needed to assemble the heatmap
#     01. GroupColorMatrix:
#     02. DifExpMatx:
#     03. colv: column dendrogram object
#     04. rowv: row dendrogram object
#     05. break_seq:
#     06. label_rows:
#     07. label_cols:
#     08. output_directory:
#     09. DistanceMethod:
#     10. ClusterMethod:
#     11. C_col: column hclust object
#     12. C_row: row hclust object
#     13. Cor_col: column correlation matrix
#     14. Cor_row: row correlation matrix
#     15. presentation:
#     16. n_clusters: number of clusters to divide the row dentrogram into
#     17. graphics_type:


{
    #Extract the data required to make a heatmap
    data_location <- dl
    #ADR: ArrangeData Return
    ADR <- ArrangeData(ColGroupsScheme,
                       replicate_scheme,
                       transformation,
                       data,
                       data_location,
                       select_rows,
                       select_groups,
                       visualization=visualization,
                       ddt,
                       handle_blanks,
                       inclusion_grouping_scheme,
                       ttest=FALSE,
                       select_rows_after_transform,
                       transform_after_column_exclusion,
                       FilterRowsForMeanSpread)
    browser()

    sig_test_list <- ADR[[1]]
    output_directory <- ADR[[2]]
    group_order <- ADR[[3]]
    gene_name <- ADR[[4]]
    DATA_long <- ADR[[5]]
    FillColors <- ADR[[6]]
    DATA <- ADR[[7]]
    GroupColorMatrix <- ADR[[8]]
    groups_corresponding <- ADR[[9]]
    DATA_original <- ADR[[11]] #the data before any transformations or exclusions were made

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
    #Make Assemble Heatmap List
    AHI <- list(GroupColorMatrix,
                DifExpMatx,
                colv,
                rowv,
                break_seq,
                label_rows,
                label_cols,
                output_directory,
                DistanceMethod,
                ClusterMethod,
                C_col,
                C_row,
                Cor_col,
                Cor_row,
                presentation,
                n_clusters,
                graphics_type)

    AHR <- assemble_heatmap(AHI)
    heat_map_colors <- AHR[[1]]
    cutree_genes <- AHR[[2]]

    #Draw and save the legends
    MakeColorKey(break_seq,heat_map_colors,output_directory)
    MakeGroupLegend(groups_corresponding,GroupColorMatrix,ColGroupsScheme,output_directory)
    #MakeClusterLegend(n_clusters,cluster_colors,output_directory)

    #Return what could be used
    return(list(ADR,AHR,AHI))
}
