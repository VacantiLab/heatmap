OrderGroups <- function(select_groups,group_concationation,groups_corresponding,GroupColorMatrix,COLOR_KEY,groups_concatonated,colors_concatonated,gene_name,DATA_long,ddt,visualizaion)
{
  #specify the order in which the groups will be plotted and ensure they map to their corresponding colors
  #If there is nothing specified in select_groups, the order is by the order of appearance in group_color_key.txt
  group_order = NULL
  FillColors = NULL
  
  groups_corresponding_for_order <- groups_corresponding
  if (!is.null(ddt)){groups_corresponding_for_order <- groups_corresponding[,colnames(groups_corresponding)==ddt,drop=FALSE]}

  if (!is.null(groups_corresponding) && !group_concationation)
  {
      groups_in_plot <- as.character(unique(groups_corresponding_for_order))
      indices_to_keep <- rownames(COLOR_KEY) %in% groups_in_plot
      group_order <- rownames(COLOR_KEY)[indices_to_keep]
      group_order <- matrix(as.character(group_order))
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
  if (is.null(ddt))
  {
      if (!is.null(groups_corresponding)){DATA_long$group <- factor(DATA_long$group,group_order)} #this sets the order of the groups to match group_order
  }
  
  if (!is.null(ddt))
  {
    if (!is.null(groups_corresponding)){DATA_long$ddt_group <- factor(DATA_long$ddt_group,group_order)} #this sets the order of the groups to match group_order
  }
  
  DATA_long$gene <- factor(DATA_long$gene,gene_name) #this sets the order of the genes to match gene_names

  #assemble variables to return
  OrderGroups_return <- list(DATA_long,FillColors,group_order)
}
