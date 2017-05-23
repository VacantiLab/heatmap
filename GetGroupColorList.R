GetGroupColorList <- function(GROUP_KEY,DATA,group_color_designations_file,ColGroupsScheme)
{
  COLOR_KEY <- read.table(file=group_color_designations_file,head=TRUE,sep='\t',stringsAsFactors=FALSE,comment.char="",check.names=FALSE) #reads the colors associated with the groups, check.names=FALSE ensures text in column names is not changed
  COLOR_KEY <- COLOR_KEY[,-1] #remove the first column which contains names of row, which is just 'color'
                              #the group names are the column names (can be from multiple grouing schemes)

  n_samples <- ncol(DATA)
  n_designations <- length(ColGroupsScheme)

  color_df <- as.data.frame(matrix(NA,nrow=n_designations,ncol=n_samples)) #initialize a list that is going to carry the colors associated with each sample
  rownames(color_df) <- ColGroupsScheme
  colnames(color_df) <- colnames(DATA)
  group_name_df <- color_df #initialize a list that is going to carry the group name associated with each sample

  for (i in 1:n_samples) #iterate though each sample in your data
  {
    sample_name <- colnames(DATA)[i]
    group_name_df[,sample_name] <- GROUP_KEY[,sample_name]
    current_groups <- group_name_df[,sample_name]
    color_df[,sample_name] <- as.character(COLOR_KEY[current_groups])
  }
  #make the requisite matrices out of the above made lists
  color_matrix <- as.matrix(color_df)
  group_name_matrix <- as.matrix(group_name_df)

  GroupColorListReturn <- list(color_matrix,group_name_matrix) #create list to return all that is required by MakeHeatmap() function
  return(GroupColorListReturn)
}
