MakeScatterLinePlot <- function(data_location,
                                ColGroupsScheme=NULL,
                                transformation=NULL,
                                data=NULL,
                                select_rows=NULL,
                                select_rows_after_transform = NULL,
                                select_groups=NULL,
                                inclusion_grouping_scheme=NULL,
                                replicate_scheme=NULL,
                                ddt=NULL,
                                med_norm=FALSE,
                                handle_blanks='remove_row',
                                transform_after_column_exclusion = FALSE,
                                plot_width = 15,
                                plot_height = 5)
#This function lists samples across the x axis and expression along the y
#Each gene is treated as a group
# Inputs:
# data_location: a pathway to where the text file containing the data is stored, must have '/' at the end
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
# select_groups: This can be an array of the group names that are to be plotted or a list of arrays of group names or NULL
#    If it is an array of group names, those groups are the only ones plotted
#    If it is a list of arrays of group names, groups in the same array are combined into a single group
#    If it is NULL, all groups in the ColGroupsScheme are plotted
#    These can be groups outside of ColGroupsScheme, but the scheme must then be specified as the inclusion_grouping_scheme
#    The spcification of transform_after_column_exclusion determines if this impacts data transformation/normalization

{
  #Extract the data required to make a scatter plot
  ArrangeData_return <- ArrangeData(ColGroupsScheme,
                                    replicate_scheme,
                                    transformation,
                                    data,
                                    data_location,
                                    select_rows,
                                    select_groups,
                                    visualization='ScatterLinePlot',
                                    ddt=NULL,
                                    handle_blanks,
                                    inclusion_grouping_scheme,
                                    ttest=FALSE,
                                    select_rows_after_transform,
                                    transform_after_column_exclusion)
  sig_test_list <- ArrangeData_return[[1]]
  output_directory <- ArrangeData_return[[2]]
  group_order <- ArrangeData_return[[3]]
  gene_name <- ArrangeData_return[[4]]
  DATA_long <- ArrangeData_return[[5]]
  FillColors <- ArrangeData_return[[6]]
  DATA <- ArrangeData_return[[7]]
  GroupColorMatrix <- ArrangeData_return[[8]]
  groups_corresponding <- ArrangeData_return[[9]]
  DATA_original <- ArrangeData_return[[11]]

  #Make the scatter plot
  b <- assemble_scatter_line_plot(DATA_long,GroupColorMatrix,output_directory,select_rows,plot_width,plot_height)

  MSLPR <- list(DATA,DATA_original)
  return(MSLPR)
}
