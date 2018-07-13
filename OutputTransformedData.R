OutputTransformedData <- function(data_location,ColGroupsScheme=NULL,transformation=NULL,replicate_scheme=NULL,data=NULL,select_rows=NULL,select_groups=NULL,ddt=NULL,med_norm=FALSE)

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
#    If it is a path to a .txt file, that file must contain the gene names, one on each line, and they will be used as the select_rows
# select_groups: This can be an array of the group names (from the provided ColGroupsScheme) that are to be plotted or a list of arrays of group names or NULL
#    If it is an array of group names, those groups are the only ones plotted
#    If it is a list of arrays of group names, groups in the same array are combined into a single group
#    The ratio plotted is the first group median over the second group median
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
# med_norm specifies to median normalize columns, will occur before rows are selected or other transformations performed (Done within OpenDataFile)

{
    #Extract the data required to make a heatmap
    ArrangeData_return <- ArrangeData(ColGroupsScheme,replicate_scheme,transformation,data,data_location,select_rows,select_groups,visualization='none',ddt,med_norm)
    sig_test_list <- ArrangeData_return[[1]]
    output_directory <- ArrangeData_return[[2]]
    group_order <- ArrangeData_return[[3]]
    gene_name <- ArrangeData_return[[4]]
    DATA_long <- ArrangeData_return[[5]]
    FillColors <- ArrangeData_return[[6]]
    DATA <- ArrangeData_return[[7]]
    GroupColorMatrix <- ArrangeData_return[[8]]
    groups_corresponding <- ArrangeData_return[[9]]

    #Return what could be used
    return(list(DATA,groups_corresponding,GroupColorMatrix))

  }
