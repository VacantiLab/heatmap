MakeVolcanoPlot <- function(data_location,ColGroupsScheme=NULL,transformation=NULL,data=NULL,select_rows=NULL,select_groups=NULL,replicate_scheme=NULL,genes_to_label=NULL,med_norm=FALSE,ddt=NULL)
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
#    or it is a string that is the path to a text file (ending in .txt) that is a column list of the genes to label
# ddt: data dependent transform
#    Allows for the transformation of the columns to ratios of columns
#    It is a list of two arrays
#        The first array contains the names of the columns of the original DATA set that will be the numerators
#        The second array contains the names of the columns of the original DATA set that will be the denominators
#    It can also be a string that is one of the members of ColGroupsScheme
#        If this is the case, the data within groups designated by the ddt grouping scheme will be normalized to the median of the group

{
    #Stop the program if the transformation is log2
    CheckStop(7,parameters=list(transformation))

    #Extract the data required to make a volcano plot
    ArrangeData_return <- ArrangeData(ColGroupsScheme,replicate_scheme,transformation,data,data_location,select_rows,select_groups,visualization='volcanoplot',ddt,med_norm)
    sig_test_list <- ArrangeData_return[[1]]
    output_directory <- ArrangeData_return[[2]]
    groups_corresponding <- ArrangeData_return[[9]]
    ColGroupsScheme <- ArrangeData_return[[11]]



    #Designate the data visualized in the volcano plot
    volcano_df <- sig_test_list[[1]] #there is only one member of this list for a volcano plot because there are must only be two grouping schemes compared

    #Make the plot
    genes_to_label <- GetGeneList(genes_to_label,data_location)
    vp <- assemble_volcano_plot(volcano_df,output_directory,genes_to_label,XData='log2_ratio',YData='nlog10_p',filename='volcano.pdf')

    #Rank the genes in order of over-expression
    RankVolcanoData_return <- RankVolcanoData(volcano_df,output_directory)
    ranked_volcano_df <- RankVolcanoData_return[[1]]
    lin_model <- RankVolcanoData_return[[2]]

    #Make the plot where the rankings were derived from
    rvp <- assemble_volcano_plot(ranked_volcano_df,output_directory,genes_to_label,XData='reg_log2_ratio',YData='reg_nlog10_p',filename='volcano_regression.pdf')

    return(sig_test_list)
}
