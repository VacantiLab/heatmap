RetrieveGroups <- function(DATA,ColGroupsScheme,group_designations_file,group_color_designations_file,select_groups,replicate_scheme)
{
    #Determine if there is column or row grouping
    ColGroups <- FALSE
    groups_corresponding <- NULL
    GroupColorMatrix <- NULL
    COLOR_KEY <- NULL
    if (is.character(ColGroupsScheme)) {ColGroups <- TRUE} #this will be true if the input of ColGroupsScheme was NULL and there is a replicate_scheme

    if (ColGroups)
    {
        #Assign group memberships (ie: those that grow, those that live, those in HER2 cluster, etc.)
        #Read the file containing the group designations
        GROUP_KEY <- read.table(file=group_designations_file,head=TRUE,check.names=FALSE,sep='\t',stringsAsFactors=FALSE) #check.names=FALSE prevents changing special characters
        rownames(GROUP_KEY) <- GROUP_KEY[,1]
        possible_group_schemes <- rownames(GROUP_KEY)
        all_ColGroupScheme_real <- CheckAllIn(ColGroupsScheme,possible_group_schemes)
        if (!all_ColGroupScheme_real){stop('custom message: a specified ColGroupsScheme or the replicate_scheme does not exist.')}
        GROUP_KEY <- GROUP_KEY[,-1] #remove the first column which contains the name of the group designation system (i.e. PAM50, Protein Clustering, etc.)
        GROUP_KEY <- GROUP_KEY[ColGroupsScheme,] #the group key rows corresponding to the grouping schemes considered are selected
        #Get the vector of colors corresponding to the group membership of each patient
        GroupColorListReturn <- GetGroupColorList(GROUP_KEY,DATA,group_color_designations_file,ColGroupsScheme)
        GroupColorMatrix <- GroupColorListReturn[[1]]
        groups_corresponding <- GroupColorListReturn[[2]]
        COLOR_KEY <- GroupColorListReturn[[3]]
    }

    #return used variables
    RetrieveGroups_return <- list(groups_corresponding,GroupColorMatrix,COLOR_KEY)
}

###############################################################################

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
  color_matrix <- t(as.matrix(color_df))
  group_name_matrix <- t(as.matrix(group_name_df))

  GroupColorListReturn <- list(color_matrix,group_name_matrix,COLOR_KEY) #create list to return all that is required by MakeHeatmap() function
  return(GroupColorListReturn)
}

###############################################################################

SelectGroups <- function(select_groups,DATA,ColGroupsScheme,groups_corresponding,GroupColorMatrix,inclusion_grouping_scheme)
#Select the specified group names (if applicable)
{
    if (is.character(select_groups))
    {
        select_groups_in_same_scheme <- CheckAllIn(select_groups,groups_corresponding[,inclusion_grouping_scheme])
        if (!select_groups_in_same_scheme){stop('custom message: Not all select_groups are in the same ColGroupsScheme. You cannot specify to consider only group members from multiple grouping schemes.')}
        #not really relevant for boxplots - check for heatmaps
        DATA2 <- DATA #make a copy of the DATA data frame because appending non-numeric rows changes the numeric class of its contents
        DATA2[ColGroupsScheme,sort(colnames(DATA))] <- t(groups_corresponding[sort(colnames(DATA)),ColGroupsScheme]) #add rows corresponding to the group classifications of each sample
        DATA2[paste(ColGroupsScheme,'- color'),sort(colnames(DATA))] <- t(GroupColorMatrix[sort(colnames(DATA)),ColGroupsScheme]) #add rows corresponding to the group color assignments of each sample
        col_to_keep_indices <- DATA2[inclusion_grouping_scheme,] %in% select_groups #find the columns to keep - this gives the indices
        DATA2 <- DATA2[,col_to_keep_indices] #keep the specified columns, this is done in DATA2 so the GroupColorMatrix and groups_corresponding matrix can be updated
        DATA <- DATA[,col_to_keep_indices] #keep the specified columns in the data frame to be plotted
        GroupColorMatrix <- t(as.matrix(DATA2[paste(ColGroupsScheme,'- color'),])) #update the GroupColorMatrix
        colnames(GroupColorMatrix) <- ColGroupsScheme #restore the rownames to their original (without the '- color')
        groups_corresponding <- t(as.matrix(DATA2[ColGroupsScheme,])) #update the groups_corresponding matrix
        rm(DATA2) #remove the now unnecessary DATA2 data frame from memory
    }

    SelectGroups_return <- list(DATA,groups_corresponding,GroupColorMatrix)
    return(SelectGroups_return)
}

###############################################################################

ConcatonateGroups <- function(group_divisions,groups_corresponding,GroupColorMatrix,COLOR_KEY)
# Inputs:
# group_divisions is a list of arrays where each array is a group of groups to be concatonated
# groups_corresponding is a column matrix containing the un-concatonated group assignments for each sample corresponding to the rows
# GroupColorMatrix is a column matrix containing the un-concatonated color assignments for each sample corresponding to the rows
# COLOR_KEY is a data frame containing the group names across the columns and in its single row are the corresponding colors
#     These are the color assignments for unique group names, NOT for samples
# Outputs:
# groups_corresponding is a column matrix containing the concatonated group assignments for each sample corresponding to the rows
# GroupColorMatrix is a column matrix containing the concatonated group assignments for each sample corresponding to the rows
# groups_concatonated is a column matrix containing the unique group names, ordered as they were input in group_divisions
# colors_concatonated is a column matrix containing the unique color assignments, ordered as corresponding to groups_concatonated
{
    group_concatonation = is.list(group_divisions)
    groups_concatonated = NULL
    colors_concatonated = NULL

    if (group_concatonation)
    {
        #get a list of the new group names
        groups_concatonated <- lapply(group_divisions,paste,collapse=':')

        #make that list an array
        n_groups_concatonated <- length(groups_concatonated)
        groups_concatonated <- matrix(unlist(groups_concatonated),nrow=n_groups_concatonated)

        #replace the group names with their corresponding concatonated names
        #do likewise for the colors
        colors_concatonated <- matrix(ncol=1,nrow=n_groups_concatonated)
        rownames(colors_concatonated) <- groups_concatonated
        for (i in 1:n_groups_concatonated)
        {
            current_group_concatonated <- groups_concatonated[i]
            concatonate_indices <- groups_corresponding %in% group_divisions[[i]]
            groups_corresponding[concatonate_indices] <- groups_concatonated[i]
            GroupColorMatrix[concatonate_indices,1] <- COLOR_KEY[1,group_divisions[[i]][1]] #The corresponding color for each contatonated group is the corresponding color to the first member sub-group
            colors_concatonated[current_group_concatonated,1] <- COLOR_KEY[1,group_divisions[[i]][1]]
        }
     }

     #return desired variables
     ConcatonateGroups_return <- list(groups_corresponding,GroupColorMatrix,groups_concatonated,colors_concatonated)
     return(ConcatonateGroups_return)
}

###############################################################################

MedianGroup <- function(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme,ColGroupsScheme)
{
#If data is presented in duplicates or triplicates, one may want to consider the medians while clustering
#Replicates must be specified as a grouping scheme, and this replicates grouping scheme
if (is.character(replicate_scheme))
{
    StatTransformByGroup_return <- StatTransformByGroup(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme,ColGroupsScheme)
    DATA <- StatTransformByGroup_return[[1]]
    groups_corresponding <- StatTransformByGroup_return[[2]]
    GroupColorMatrix <- StatTransformByGroup_return[[3]]
}

#If the replicate scheme was the only ColGroupsScheme, then the returned groups_corresponding and GroupColorMatrix will be empty column matrices
if (length(groups_corresponding)==0)
    {
        groups_corresponding = NULL
        GroupColorMatrix = NULL
    }

MedianGroup_return <- list(DATA,groups_corresponding,GroupColorMatrix)
return(MedianGroup_return)
}

###############################################################################

StatTransformByGroup <- function(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme,ColGroupsScheme)
#This function finds the median of the specified group in replicate scheme and creates a new data frame with the median values
{
    groups_corresponding_rep_scheme <- groups_corresponding[,replicate_scheme,drop=FALSE]
    unique_groups_rep_scheme <- unique(groups_corresponding_rep_scheme)
    n_unique_groups_rep_scheme <- length(unique_groups_rep_scheme)
    samples <- colnames(DATA)

    DATA2 <- data.frame(matrix(ncol=0,nrow=nrow(DATA)))
    rownames(DATA2) <- rownames(DATA)

    ColGroupsScheme <- ColGroupsScheme[ColGroupsScheme!=replicate_scheme]
    GroupColorMatrix <- GroupColorMatrix[,!colnames(GroupColorMatrix)==replicate_scheme,drop=FALSE]
    groups_corresponding <- groups_corresponding[,!colnames(groups_corresponding)==replicate_scheme,drop=FALSE]

    n_ColGroupsScheme <- ncol(groups_corresponding) #these are the remaining grouping schemes (other than that used to denote replicates)

    Color_DF <- data.frame(matrix(nrow=n_unique_groups_rep_scheme,ncol=n_ColGroupsScheme))
    rownames(Color_DF) <- unique_groups_rep_scheme #these are the new sample names (the replicates are grouped to samples)
    colnames(Color_DF) <- ColGroupsScheme #these are the remaining grouping schemes of the grouped-replicate samples
    Groupings_DF <- data.frame(matrix(nrow=n_unique_groups_rep_scheme,ncol=n_ColGroupsScheme))
    rownames(Groupings_DF) <- unique_groups_rep_scheme #these are the new sample names (the replicates are grouped to samples)
    colnames(Groupings_DF) <- ColGroupsScheme #these are the remaining grouping schemes of the grouped-replicate samples

    for (i in 1:n_unique_groups_rep_scheme)
    {
        group_rep <- unique_groups_rep_scheme[i]
        member_indices <- groups_corresponding_rep_scheme==group_rep
        members <- rownames(groups_corresponding_rep_scheme[member_indices,,drop=FALSE])
        median_vector <- apply(DATA[,members],1,median)
        DATA2[,group_rep] <- median_vector
        group_color_row <- which(member_indices)[1]
        Color_DF[group_rep,] <- GroupColorMatrix[group_color_row,]
        Groupings_DF[group_rep,] <- groups_corresponding[group_color_row,]
    }

    GroupColorMatrix <- as.matrix(Color_DF)
    groups_corresponding <- as.matrix(Groupings_DF)

    DATA <- DATA2

    StatTransformByGroup_return <- list(DATA,groups_corresponding,GroupColorMatrix)
    return(StatTransformByGroup_return)
}
