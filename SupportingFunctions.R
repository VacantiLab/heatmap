OpenDataFile <- function(data,select_rows)
#Takes a data file location and reads the data into a data frame compatible with the MakeHeatMap function
{
    #Import the data
    if (is.character(data)) {DATA <- read_txt_to_df(data)} #If the input, data, is provided as a string it is the directory to the text file with the data
    if (is.data.frame(data)) {DATA <- data} #If the input, data, is provided as a data frame, it is the data

    #remove any rows with NA as an entry
    has_no_na_row_indices <- apply(DATA,1,NoNA)
    DATA <- DATA[has_no_na_row_indices,]

    #Select rows from the data frame
    if (!is.null(select_rows)) {DATA <- DATA[rownames(DATA) %in% select_rows,]}
    #selecting the rows by name and not position dictates the order of the rows, thus the rows are selected by position here to maintain the order which
    #is important in consitently arranging equivalent positions in the dendrogram below. Equivalent positions are two members linked at the lowest possible level.

    return(DATA)
}

###############################################################################

StoreHeatmap <- function()
{
    #Input where the heatmap is stored
    working_directory <- getwd()
    working_directory_up1 <- gsub('/[^/]*$','/',working_directory) #matches '/' followed by 0 or more characters other than '/' followed by the end of the string, and replaces with '/'
    HeatmapDirectory <- paste(working_directory_up1,'output/',sep='')
    dir.create(HeatmapDirectory) #creates the directory for the output
    return(HeatmapDirectory)
}

###############################################################################

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
        inclusion_grouping_scheme_indices <- apply(DATA2[ColGroupsScheme,], 1, function(r) any(r %in% select_groups_to_keep)) #find which group classification scheme is being used to make the selection - this gives the indices of the rows in DATA2
        inclusion_grouping_scheme <- rownames(DATA2[ColGroupsScheme,])[inclusion_grouping_scheme_indices] #this gives the classification scheme used to make the selection of groups to include
        col_to_keep_indices <- DATA2[inclusion_grouping_scheme,] %in% select_groups_to_keep #find the columns to keep - this gives the indices
        DATA2 <- DATA2[,col_to_keep_indices] #keep the specified columns, this is done in DATA2 so the GroupColorMatrix and groups_corresponding matrix can be updated
        DATA <- DATA[,col_to_keep_indices] #keep the specified columns in the data frame to be plotted
        GroupColorMatrix <- as.matrix(DATA2[paste(ColGroupsScheme,'- color'),]) #update the GroupColorMatrix
        rownames(GroupColorMatrix) <- ColGroupsScheme #restore the rownames to their original (without the '- color')
        groups_corresponding <- as.matrix(DATA2[ColGroupsScheme,]) #update the groups_corresponding matrix
        rm(DATA2) #remove the now unnecessary DATA2 data frame from memory
    }

    #If data is presented in duplicates or triplicates, one may want to consider the medians while clustering
    #Replicates must be specified as a grouping scheme, and this replicates grouping scheme
    if (is.character(replicate_scheme))
    {
        StatTransformByGroup_return <- StatTransformByGroup(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme)
        DATA <- StatTransformByGroup_return[[1]]
        groups_corresponding <- StatTransformByGroup_return[[2]]
        GroupColorMatrix <- StatTransformByGroup_return[[3]]
    }

    #The group color matrix for the columns is returned with columns corresponding to columns, this is the transpose of what it needs to be
    GroupColorMatrix <- t(GroupColorMatrix)

    #return used variables
    RetrieveGroups_return <- list(DATA,groups_corresponding,GroupColorMatrix,COLOR_KEY)
}

###############################################################################

ClusterData <- function(DATA,DistanceMethod,ClusterMethod,rev_c_dend)
{
    #compute the distance matrix depending on the desired method (euclidan, pearson, or spearman)
    if (DistanceMethod == 'euclidian')
    {
        D_col <- dist(t(DATA)) #This function computes the distance matrix where the rows are points, but in this instance the columns should be points, so use the transpose of the data matrix.
        D_row <- dist(DATA) #This function computes the distance matrix where the rows are points. This is the desired output in this instance.
    }

    if (DistanceMethod == 'pearson' || DistanceMethod == 'spearman')
    {
        Cor_col <- cor(DATA,method=DistanceMethod) #transpose oppsitely as when using dist() function?
        D_col <- (1-Cor_col)/2 #this is actually dissimilarity, zero is least dissimilar (correlation of 1) and 1 is most dissimilar (correlation of -1)
        D_col <- as.dist(D_col) #make a distance object so it can be input into hclust, clustering is performed on dissimilarity
        Cor_row <- cor(t(DATA),method=DistanceMethod)
        D_row <- (1-Cor_row)/2
        D_row <- as.dist(D_row)
    }
    #This function computes the distance matrix where the rows are points. This is the desired output in this instance.

    #A distance matrix looks like the following where 1->2 indicates distance from the point
    #indicated by row 1 to that indicated by row 2. It is a lower left triangular matrix. A
    #distance matrix is not actually a matrix object, rather a distance object.
    # 1->2  ---- ---- ----
    # 1->3  2->3 ---- ----
    # 1->4  2->4 3->4 ----
    # 1->5  2->5 3->5 4->5
    #Note the function dist() is not to be confused with as.dist(). as.dist() takes a square matrix
    #and returns the lower triangular portion as a distance object.

    #A correlation matrix looks like the following where 1->2 indicates distance from the point
    #indicated by column 1 to that indicated by column 2. It is a lower left triangular matrix.
    #Use as.dist to convert a correlation matrix to a distance object by taking the lower triangular matrix.
    # 1->1  2->1  3->1  4->1
    # 1->2  2->2  3->2  4->2
    # 1->3  2->3  3->3  4->3
    # 1->4  2->4  3->4  4->4
    # 1->5  2->5  3->5  4->5

    #cluster the data by feeding the distance matrices to hclust().
    #This clusters the columns of the data matrix
    C_col <- hclust(D_col, method=ClusterMethod,members=NULL)
    #This clusters the rows of the data matrix
    C_row <- hclust(D_row, method=ClusterMethod,members=NULL)

    #create a dendrogram objects to specify the order and dendrogram labeling of the columns and rows in the heatmap
    colv <- as.dendrogram(C_col)
    rowv <- as.dendrogram(C_row)
    #the default places the higher "next level" branch point on the right side of a branch point. At the end of the tree where the end points are
    #of equal heights, and the right/left position is dictated by position in the dist object (ultimately dictated by position in the original data frame)

    #reverse colum dendrogram if specified to do so
    if (rev_c_dend) {colv <- rev(colv)}

    #return used variables
    ClusterData_return <- list(colv,rowv,C_col,C_row)
}

###############################################################################

NoNA <- function(vector)
#Function to return TRUE if the row does not have any NAs
{
    NoNA <- !is.na(sum(vector))
}

###############################################################################

read_txt_to_df <- function(txt_directory)
#Function to read a tab delimited text file into a data frame
{
    DF <- read.table(file=txt_directory,head=TRUE,check.names=FALSE,sep='\t') #check.names=FALSE prevents an 'X' from being added to the numeric column names
    #Name the rows of the data frame as the genes given in the first column of the data frame
    RowNames <- as.character(DF[,1])
    rownames(DF) <- RowNames
    DF[,1] <- NULL #remove the first column of the data frame as it is no longer needed
    return(DF)
}

###############################################################################

FindRowLabels <- function(label_rows,DifExpMatx)
#Function to return row label_rows
{
    #Provide all of the rownames if label_rows input is TRUE
    if (is.logical(label_rows))
    {
        if(label_rows){label_rows <- rownames(DifExpMatx)}
    }

    #If there are rownames specified to be labeled, ensure those are the only ones that are labeled
    if (is.character(label_rows))
    {
        heatmap_rownames <- rownames(DifExpMatx)
        heatmap_rowname_indices_to_remove <- !(heatmap_rownames %in% label_rows)
        heatmap_rownames[heatmap_rowname_indices_to_remove] <- ''
        label_rows <- heatmap_rownames
    }
return(label_rows)
}

###############################################################################

FindColLabels <- function(label_cols,DifExpMatx)
#Function to return label_cols
{
    #Provide all of the colnames if label_cols input is TRUE
    if (is.logical(label_cols))
    {
        if(label_cols){label_cols <- colnames(DifExpMatx)}
    }

    #If there are colnames specified to be labeled, ensure those are the only ones that are labeled
    if (is.character(label_cols))
    {
        heatmap_colnames <- colnames(DifExpMatx)
        heatmap_colname_indices_to_remove <- !(heatmap_colnames %in% label_cols)
        heatmap_colnames[heatmap_colname_indices_to_remove] <- ''
        label_cols <- heatmap_colnames
    }
return(label_cols)
}
