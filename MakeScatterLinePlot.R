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

{
  #Extract the data required to make a scatter plot
  ArrangeData_return <- ArrangeData(ColGroupsScheme,replicate_scheme,transformation,data,data_location,select_rows,select_groups,visualization='ScatterLinePlot',ddt=NULL,handle_blanks,inclusion_grouping_scheme,ttest=FALSE,select_rows_after_transform,transform_after_column_exclusion)
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
