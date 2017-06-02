RetrieveGroups <- function(DATA,ColGroupsScheme,group_designations_file,group_color_designations_file,select_groups,replicate_scheme)
{
    #Determine if there is column or row grouping
    ColGroups <- FALSE
    if (is.character(ColGroupsScheme)) {ColGroups <- TRUE}

    if (ColGroups || RowGroups)
    {
        #Assign group memberships (ie: those that grow, those that live, those in HER2 cluster, etc.)
        #Read the file containing the group designations
        GROUP_KEY <- read.table(file=group_designations_file,head=TRUE,check.names=FALSE,sep='\t',stringsAsFactors=FALSE) #check.names=FALSE prevents changing special characters
        rownames(GROUP_KEY) <- GROUP_KEY[,1]
        GROUP_KEY <- GROUP_KEY[,-1] #remove the first column which contains the name of the group designation system (i.e. PAM50, Protein Clustering, etc.)
        GROUP_KEY <- GROUP_KEY[ColGroupsScheme,] #the group key rows corresponding to the grouping schemes considered are selected
        #Get the vector of colors corresponding to the group membership of each patient
        GroupColorListReturn <- GetGroupColorList(GROUP_KEY,DATA,group_color_designations_file,ColGroupsScheme)
        GroupColorMatrix <- GroupColorListReturn[[1]]
        groups_corresponding <- GroupColorListReturn[[2]]
        COLOR_KEY <- GroupColorListReturn[[3]]
        #group_colors <- GroupColorListReturn[[3]]
        #group_names <- GroupColorListReturn[[4]]
    }

    #Select the specified group names (if applicable)
    #This must occur before scaling below!
    if (is.character(select_groups))
    {
        DATA2 <- DATA #make a copy of the DATA data frame because appending non-numeric rows changes the numeric class of its contents
        DATA2[ColGroupsScheme,sort(colnames(DATA))] <- groups_corresponding[ColGroupsScheme,sort(colnames(DATA))] #add rows corresponding to the group classifications of each sample
        DATA2[paste(ColGroupsScheme,'- color'),sort(colnames(DATA))] <- GroupColorMatrix[ColGroupsScheme,sort(colnames(DATA))] #add rows corresponding to the group color assignments of each sample
        inclusion_grouping_scheme_indices <- apply(DATA2[ColGroupsScheme,], 1, function(r) any(r %in% select_groups)) #find which group classification scheme is being used to make the selection - this gives the indices of the rows in DATA2
        inclusion_grouping_scheme <- rownames(DATA2[ColGroupsScheme,])[inclusion_grouping_scheme_indices] #this gives the classification scheme used to make the selection of groups to include
        col_to_keep_indices <- DATA2[inclusion_grouping_scheme,] %in% select_groups #find the columns to keep - this gives the indices
        DATA2 <- DATA2[,col_to_keep_indices] #keep the specified columns, this is done in DATA2 so the GroupColorMatrix and groups_corresponding matrix can be updated
        DATA <- DATA[,col_to_keep_indices] #keep the specified columns in the data frame to be plotted
        GroupColorMatrix <- as.matrix(DATA2[paste(ColGroupsScheme,'- color'),]) #update the GroupColorMatrix
        rownames(GroupColorMatrix) <- ColGroupsScheme #restore the rownames to their original (without the '- color')
        groups_corresponding <- as.matrix(DATA2[ColGroupsScheme,]) #update the groups_corresponding matrix
        rm(DATA2) #remove the now unnecessary DATA2 data frame from memory
    }

    #If data is presented in duplicates or triplicates, one may want to consider the medians while clustering
    #Replicates must be specified as a grouping scheme, and this replicates grouping scheme
    remaining_grouping_scheme <- TRUE #whether there are grouping schemes besides one used to mark replicates
    if (is.character(replicate_scheme))
    {
        StatTransformByGroup_return <- StatTransformByGroup(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme)
        DATA <- StatTransformByGroup_return[[1]]
        groups_corresponding <- StatTransformByGroup_return[[2]]
        GroupColorMatrix <- StatTransformByGroup_return[[3]]
        group_values <- StatTransformByGroup_return[[4]]
        remaining_grouping_scheme <- StatTransformByGroup_return[[5]]
    }

    #The group color matrix for the columns is returned with columns corresponding to columns, this is the transpose of what it needs to be
    if (remaining_grouping_scheme){GroupColorMatrix <- t(GroupColorMatrix)}

    #return used variables
    RetrieveGroups_return <- list(DATA,groups_corresponding,GroupColorMatrix,COLOR_KEY,group_values)
}

###############################################################################

ConcatonateGroups <- function(group_divisions,groups_corresponding,GroupColorMatrix,COLOR_KEY)
{
    group_concatonation = is.list(group_divisions)

    if (group_concatonation)
    {
        #get a list of the new group names
        groups_concatonated <- lapply(group_divisions,paste,collapse=':')

        #make that list an array
        groups_concatonated <- unlist(groups_concatonated)

        #replace the group names with their corresponding concatonated names
        #do likewise for the colors
        n_groups_concatonated <- length(groups_concatonated)
        for (i in 1:n_groups_concatonated)
        {
            concatonate_indices <- groups_corresponding %in% group_divisions[[i]]
            groups_corresponding[concatonate_indices] <- groups_concatonated[i]
            GroupColorMatrix[concatonate_indices,1] <- COLOR_KEY[1,group_divisions[[i]][1]] #The corresponding color for each contatonated group is the corresponding color to the first member sub-group
        }
     }

     #return desired variables
     ConcatonateGroups_return <- list(groups_corresponding,GroupColorMatrix)
     return(ConcatonateGroups_return)
}

###############################################################################
