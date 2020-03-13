MakeBoxPlot <- function(data_location,ColGroupsScheme=NULL,box_plot_type='bar_plot',plot_width=2.0,plot_height=1.5,bar_width=0.75,legend_position=c(0.8,0.92),text_angle=0,ttest=FALSE,transformation=NULL,data=NULL,select_rows=NULL,select_groups=NULL,inclusion_grouping_scheme=NULL,replicate_scheme=NULL,qc_plot=FALSE,med_norm=FALSE,ddt=NULL,handle_blanks='remove_row',violin=FALSE,output_in_dl=TRUE)
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
#    If it is NULL, all genes are plotted
#    If it is a path to a .txt file, that file must contain the gene names, one on each line, and they will be used as the select_rows
# select_groups: This can be an array of the group names (from the provided ColGroupsScheme) that are to be plotted or a list of arrays of group names or NULL
#    If it is an array of group names, those groups are the only ones plotted
#    If it is a list of arrays of group names, groups in the same array are combined into a single group
#    If it is NULL, all groups in the ColGroupsScheme are plotted
# inclusion_grouping_scheme: This is the grouping scheme that you want to specify to select the columns with
#    It can remain as the default NULL and the columns will be selected based on the first grouping scheme in ColGroupsScheme
# replicate_scheme: This specifies the grouping scheme that is used to specify groups of replicates
#    This must NOT be a member of ColGroupsScheme, though it must be a grouping scheme defined in group_key.txt
#        As such each member of this grouping scheme must also have colors specified in group_color_key.txt
#    If this is specified, all members of a single group are treated as a single sample and the median values are used
# box_plot_type: This specifies the representation of box-plot like data. Inputs can be 'boxplot', 'scatter_bar_plot', 'bar_plot'.
# plot_width and plot_height are in inches
# bar_width is in units on the x-axis. For discrete units, the distance between units is considered 1. bar_width specifies the combined width of the bars at a discrete unit.
# the legend position is specified in relative cartesian coordinates in the plot area
# ttest dictates whether to perform pairwise ttests. Will give an error if a select_genes has a gene not in the dataset

{
    #Extract the data required to make a box plot
    ArrangeData_return <- ArrangeData(ColGroupsScheme,replicate_scheme,transformation,data,data_location,select_rows,select_groups,visualization='boxplot',ddt,med_norm,handle_blanks,inclusion_grouping_scheme,ttest)
    sig_test_list <- ArrangeData_return[[1]]
    output_directory <- ArrangeData_return[[2]]
    group_order <- ArrangeData_return[[3]]
    gene_name <- ArrangeData_return[[4]]
    DATA_long <- ArrangeData_return[[5]]
    FillColors <- ArrangeData_return[[6]]
    DATA_transformed <- ArrangeData_return[[7]]
    DATA_transformed_full <- ArrangeData_return[[10]]

    #Find the y-limits for the boxplot based on the data
    y_bounds <- get_y_bounds(group_order,gene_name,DATA_long)

    if (output_in_dl){output_directory <- dl}

    #Make the plot
    if (!violin){b <- assemble_box_plot(DATA_long,FillColors,output_directory,y_bounds,qc_plot,box_plot_type,plot_width,plot_height,bar_width,legend_position,text_angle)}
    if (violin){b <- assemble_violin_plot(DATA_long,FillColors,output_directory,y_bounds,qc_plot)}

    #assemble variables to return
    MakeBoxPlot_return <- list(sig_test_list,DATA_transformed,DATA_transformed_full,DATA_long)

    return(MakeBoxPlot_return)
}
