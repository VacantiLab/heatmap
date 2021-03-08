MakeScatterPlot <- function(data_location,ColGroupsScheme=NULL,transformation=NULL,data=NULL,select_rows=NULL,select_groups=NULL,replicate_scheme=NULL,ddt=NULL,med_norm=FALSE,handle_blanks='remove_row')
# select_rows_after_transform are the rows that are selected after the transformation
#     can be null for no rows selected after the transformation
{
  #Extract the data required to make a scatter plot
  ArrangeData_return <- ArrangeData(ColGroupsScheme = ColGroupsScheme,
                                    replicate_scheme = replicate_scheme,
                                    transformation = transformation,
                                    data = data,
                                    data_location = data_location,
                                    select_rows = select_rows,
                                    select_groups = select_groups,
                                    visualization='scatterplot',
                                    ddt = ddt,
                                    handle_blanks = handle_blanks,
                                    inclusion_grouping_scheme = NULL,
                                    transform_after_column_exclusion = FALSE,
                                    select_rows_after_transform = NULL,
                                    ttest = FALSE)
  
  sig_test_list <- ArrangeData_return[[1]]
  output_directory <- ArrangeData_return[[2]]
  group_order <- ArrangeData_return[[3]]
  gene_name <- ArrangeData_return[[4]]
  DATA_long <- ArrangeData_return[[5]]
  FillColors <- ArrangeData_return[[6]]
  DATA <- ArrangeData_return[[7]]
  GroupColorMatrix <- ArrangeData_return[[8]]
  groups_corresponding <- ArrangeData_return[[9]]
  
  #Make the scatter plot
  b <- assemble_scatter_plot(DATA,GroupColorMatrix,output_directory,select_rows)
}
