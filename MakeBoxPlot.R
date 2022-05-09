MakeBoxPlot <- function(data_location,
                        ColGroupsScheme=NULL,
                        box_plot_type='bar_plot',
                        plot_width=2.0,
                        plot_height=1.5,
                        bar_width=0.75,
                        legend_position=c(0.8,0.92),
                        text_angle=0,
                        ttest=FALSE,
                        transformation=NULL,
                        data=NULL,
                        select_rows=NULL,
                        select_groups=NULL,
                        inclusion_grouping_scheme=NULL,
                        replicate_scheme=NULL,
                        qc_plot=FALSE,
                        med_norm=FALSE,
                        ddt=NULL,
                        handle_blanks='remove_row',
                        violin=FALSE,
                        output_in_dl=FALSE,
                        select_rows_after_transform=NULL,
                        transform_after_column_exclusion=FALSE,
                        ybounds=NULL,
                        ytick=NULL,
                        zscore_rows=FALSE,
                        ErrorBarSize=0.75,
                        PointSize=3,
                        ErrorFile=NULL,
                        x_var = 'gene',
                        y_var = 'value',
                        color_var = 'group',
                        ratio_scheme_groups=NULL,
                        data_bounds=NULL,
                        custom_y_bounds=NULL)
# data_location: a pathway to where the text file containing the data is stored, must have '/' at the end
#    The data file must be named quantities.txt with the genes down the rows and sample names across the columns
#    There must also be a group_key.txt file with the sample names down the rows and the grouping schemes across the columns
#        The group names within each grouping scheme must not match any of those of another grouping scheme
#        Each column must have a header
#    There must also be a group_color_key.txt file with the group names down the rows and the colors beside them
#        All of the group names are listed down the rows with no indication of their scheme membership (thus the names must be unique)
#        Each column must have a header
#    These files must all be tab delimited
# ColGroupsScheme: the name(s) of the grouping scheme used, indicates which row to take from group_key.txt
#    There can only be one for the MakeBoxPlot() function
# box_plot_type
#     boxplot
#     bar_plot
#     scatter_bar_plot
#     line_plot
#         The rownames are expected to be numeric values

# transformation: This specifies how the rows should be transformed.
#    Options include: 'log2', 'median_center_iqr_norm', and 'median_norm_log2_transform', and others
# data: Is a data frame containing the data to be plotted if the data is passed as a data frame instead of through the quantities.txt file in data_location
#    This option could use more testing
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
#    For any type of boxplot
#        the ddt becomes the column grouping scheme and the non-control vs. control ratio is presented
#    Need to see at what point columns are excluded by the inclusion grouping scheme - should be BEFORE DDT when transform_after_column_exclusion is TRUE (it is before DDT)
#    Need to be able to specify the type of transformation for DDT
#        Is hard-coded in med_norm_within_groups within DataFramePar.R
# ratio_scheme_groups: Is a list that works with ddt or is NULL
#     The first member of the list specifies a column grouping scheme whose members describes how to take the ddt ratio
#     The second member of the list is the group of the grouping scheme (the 1st list member) whose mean is in the denominator of the ddt ratio 
#         These ratios are performed within each ddt grouping scheme group
# replicate_scheme: This specifies the grouping scheme that is used to specify groups of replicates
#    This must NOT be a member of ColGroupsScheme, though it must be a grouping scheme defined in group_key.txt
#        As such each member of this grouping scheme must also have colors specified in group_color_key.txt
#    If this is specified, all members of a single group are treated as a single sample and the median values are used
# box_plot_type: This specifies the representation of box-plot like data. Inputs can be 'boxplot', 'scatter_bar_plot', 'bar_plot'.
# plot_width and plot_height are in inches
# bar_width is in units on the x-axis. For discrete units, the distance between units is considered 1. bar_width specifies the combined width of the bars at a discrete unit.
# the legend position is specified in relative cartesian coordinates in the plot area
# ttest dictates whether to perform pairwise ttests. Will give an error if a select_genes has a gene not in the dataset
# ybounds: sets the lower and upper y-axis limits using an array, c(LowerLimit, UpperLimit)
#          if NULL, these values will be calculated automatically
# ytick: sets the division between tick marks on the y-axis. For example it could be 0.1, 1, or whatever other value
#        if NULLm these values will be calculated automatically
# zscore_rows: if TRUE, will tranform the rows to z scores AFTER all other transformations and column exclusions
# ErrorFile is a path to a tab delimited file with error bar magnitudes for each gene for each group
#     This is only specified if the error bar values are entered separately and not calculated from the data
#     The format of this file is as follows:
#     gene   GroupName
#     Gene1  value
#     Gene2  value
# ErrorBarSize: the thickness of the lines in the error bars (not the width of the caps)
# data_bounds: a vector of the upper and lower limit of the data. Data outside of this range will be replaced with the limits
# custom_y_bounds: NULL, or a list with 3 entries
#     The first entry is the lower limit of the y-axis
#     The second entry is the upper limit of the y-axis
#     The third is the spacing for tick marks between the lower and upper limits
{
    #Extract the data required to make a box plot
    ArrangeData_return <- ArrangeData(ColGroupsScheme = ColGroupsScheme,
                                      replicate_scheme = replicate_scheme,
                                      transformation = transformation,
                                      data = data,
                                      data_location = data_location,
                                      select_rows = select_rows,
                                      select_groups = select_groups,
                                      visualization = 'boxplot',
                                      ddt = ddt,
                                      handle_blanks = handle_blanks,
                                      inclusion_grouping_scheme = inclusion_grouping_scheme,
                                      ttest = ttest,
                                      select_rows_after_transform = select_rows_after_transform,
                                      transform_after_column_exclusion = transform_after_column_exclusion,
                                      ratio_scheme_groups = ratio_scheme_groups)
    sig_test_list <- ArrangeData_return[[1]]
    output_directory <- ArrangeData_return[[2]]
    group_order <- ArrangeData_return[[3]]
    gene_name <- ArrangeData_return[[4]]
    DATA_long <- ArrangeData_return[[5]]
    FillColors <- ArrangeData_return[[6]]
    DATA_transformed <- ArrangeData_return[[7]]
    DATA_transformed_full <- ArrangeData_return[[10]]
    DATA_original <- ArrangeData_return[[11]]
  

    #Find the y-limits for the boxplot based on the data
    if (is.null(ybounds)){ybounds <- get_y_bounds(group_order,gene_name,DATA_long)}

    if (output_in_dl){output_directory <- data_location}

    #Make the plot
    if (!violin){b <- assemble_box_plot(DATA_long,
                                        FillColors,
                                        output_directory,
                                        ybounds,
                                        qc_plot,
                                        box_plot_type,
                                        plot_width,
                                        plot_height,
                                        bar_width,
                                        legend_position,
                                        text_angle,
                                        transformation,
                                        ytick,
                                        ErrorBarSize,
                                        PointSize,
                                        ErrorFile,
                                        x_var,
                                        y_var,
                                        color_var,
                                        ddt,
                                        data_bounds=data_bounds,
                                        custom_y_bounds=custom_y_bounds)}
    if (violin){b <- assemble_violin_plot(DATA_long,FillColors,output_directory,y_bounds,qc_plot)}

    #assemble variables to return
    MakeBoxPlot_return <- list(sig_test_list,DATA_transformed,DATA_transformed_full,DATA_long,DATA_original)

    return(MakeBoxPlot_return)
}
