OrderGroups <- function(select_groups,group_concationation,groups_corresponding,GroupColorMatrix,COLOR_KEY,groups_concatonated,colors_concatonated,gene_name,DATA_long)
{
  #specify the order in which the groups will be plotted and ensure they map to their corresponding colors
  #If there is nothing specified in select_groups, the order is by the order of appearance in the data
  group_order = NULL
  FillColors = NULL
  if (!is.null(groups_corresponding) && !group_concationation)
  {
      group_order <- matrix(as.character(unique(groups_corresponding)),ncol=1)
      FillColors <- matrix(as.character(COLOR_KEY[group_order,1]),ncol=1)
  }

  #If groups to be plotted are specified, their order specifies the order they will be plotted in
  #This allows the user to control the order groups are plotted
  #If there are no group concatonations and select_groups is specified as an input
  if (is.character(select_groups) && !group_concationation)
  {
      group_order <- matrix(as.character(select_groups),ncol=1)
      FillColors <- matrix(as.character(COLOR_KEY[select_groups,1]),ncol=1)
  }

  #If there are group concatonations (the input is a list)
  if (group_concationation)
  {
      group_order <- matrix(as.character(groups_concatonated),ncol=1)
      FillColors <- matrix(as.character(colors_concatonated),ncol=1)
  }

  #Implement the group ordering
  if (!is.null(groups_corresponding)){DATA_long$group <- factor(DATA_long$group,group_order)} #this sets the order of the groups to match group_order
  DATA_long$gene <- factor(DATA_long$gene,gene_name) #this sets the order of the genes to match gene_names

  #assemble variables to return
  OrderGroups_return <- list(DATA_long,FillColors,group_order)
}
