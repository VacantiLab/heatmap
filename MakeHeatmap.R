MakeHeatMap <- function(data_location,ColGroupsScheme=FALSE,transformation=NULL,break_seq=NULL,replicate_scheme=NULL,DistanceMethod='pearson',ClusterMethod='ward.D2',RowGroupsScheme=FALSE,select_rows=NULL,select_groups_to_keep=NULL,rev_c_dend=FALSE,label_rows=TRUE,label_cols=TRUE)
#data is a data frame or the filename (including pathway) of where the data is stored
  #if data is a file name, it should be arranged in a tab delimited table such that table headings are across the top with row names in the 1st column:
    # Genes   Patient1    Patient2    Patient3    Patient4
    # Gene1   value       value       value       value
    # Gene2   value       value       value       value
  #if it is a data frame, it should have the patients as column names and the genes as row names
#group_designations_file is a file (including pathway) containing the group assignments of each sample
  #The first row of the file will be the sample (patient) names
  #The subsequent rows will be the corresponding group names.
  #Each row contains a group designation scheme, e.g. the second row is the PAM50 classification and the third the proteome-based classification
  #The group assignment names must be identical to those in the group_color_designations_file
#HeatmapDirectory is the directory where the heatmaps created by this function should be stored.
  #this includes the last '/'
#group_color_designations_file is a .txt file containing group names in the first row and the associated colors in the second row.
  #Right now it is for column grouping
  #The first column are the row lables and is removed by this function
  #The group assignment names must be identical to those in the group_designations_file file
  #There is one file for all of the group color designations. Group names from the different schemes can be included.
    #There will be problems if there are two identical group names from different schemes.
#DistanceMethod provides the distance used in the clustering. It can be 'pearson', 'spearman', or 'euclidian'
#break_seq specifies the breaks in the color-coding for the heatmap. It is a sequence provided by the seq(where_you_start,where_you_end,increment) function.
  #It is a required input.
#ColGroupsScheme specifies the row names in group_designations_file used as a column group classification scheme, if it is not of class "character" there is no column grouping
  #more than one classification scheme can be specified. The first column of the text file contains the column names.
  #The order the schemes appear are the same as the order they are listed in this vector, with the first appearing closest to the heatmap.
#replicate_scheme specifies if one of the ColGroupsScheme variables is a replicate scheme
  #if so, the median of the replicates will be taken and treated as a single sample in clustering
  #right now you cannot have a replicate scheme as the only grouping scheme
#RowGroups will be similar ColGroupsScheme, but is not implemented yet. Keep as FALSE
#transformation specifies whether and how the rows or the data will be transformed.
  #Currently the three operable options are 'median_center_iqr_norm', 'median_norm_log2_transform', and 'log2'.
  #The first centers the median of each row at 0 and divides by the row's median, the second divides each row by its median and takes the log2 of the resulting ratio.
  #The third option specifies to take the log2 of everything, where there is a cut-off in how negative it can go (-3 or whatever it's set to)
  #The default is FALSE where no scaling occurs.
#select_rows are the rows of the data frame you want considered
  #the default is NULL and the function considers all of the rows
#select_groups_to_keep is a list of group designations that you want to consider in the heatmap (corresponding to the column names)
  #if NULL, all columns are considered
  #This is not operable yet
#rev_c_dend reverses the column dendrogram if TRUE
#label_rows specifies whether the rows should be labeled. If TRUE, all rows are labeled, if FALSE, no rows are labeled.
  #If it is a vector of strings, only rows with names within the vector of strings are labeled.
#This function also prints the heatmap as a pdf into the specified folder
#This function returns a list (dendrograms - or any name you specify) with the row and column dendrograms.
#The row dendrogram is called by dendrograms[[1]]

{
library(gplots)
library(heatmap.plus)

#Input where the heatmap is stored
working_directory <- getwd()
working_directory_up1 <- gsub('/[^/]*$','/',working_directory) #matches '/' followed by 0 or more characters other than '/' followed by the end of the string, and replaces with '/'
HeatmapDirectory <- paste(working_directory_up1,'output/',sep='')
dir.create(HeatmapDirectory) #creates the directory for the output

#Input data
data <- paste(data_location,'quantities.txt',sep='')
group_designations_file <- paste(data_location,'group_key.txt',sep='')
group_color_designations_file <- paste(data_location,'group_color_key.txt',sep='')

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

#Determine if there is column or row grouping
ColGroups <- FALSE
if (is.character(ColGroupsScheme)) {ColGroups <- TRUE}
RowGroups <- FALSE
if (is.numeric(RowGroupsScheme)) {RowGroups <- TRUE}

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
  #group_colors <- GroupColorListReturn[[3]]
  #group_names <- GroupColorListReturn[[4]]
}

#Select the specified group names (if applicable)
#This must occur before scaling below!
if (is.character(select_groups_to_keep))
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

#Transform the data as specified
DATA <- transform_data(DATA,transformation)

#The rows with zero variance must be removed because likely because this causes NAs to be generated when calculating correlations
#Most of the time this won't be an issue because all samples are not likely to have identical values for a measurement
  #However in measruing 10k protein abundances across few samples, this could happen by chance
DATA <- remove_zero_var_rows(DATA)

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

#make numeric matrix out of data frame
DifExpMatx=matrix(as.numeric(unlist(DATA)),nrow=nrow(DATA))
rownames(DifExpMatx)<-rownames(DATA) #name the rows of the numeric matrix so they are transfered to the dendrogram in the heat map
colnames(DifExpMatx)<-colnames(DATA) #name the columns of the numeric matrix so they are transfered to the dendrogram in the heat map

#The group color matrix for the columns is returned with columns corresponding to columns, this is the transpose of what it needs to be
GroupColorMatrix <- t(GroupColorMatrix)

#Get the row names and column names for the heatmap
label_rows <- FindRowLabels(label_rows,DifExpMatx)
label_cols <- FindColLabels(label_cols,DifExpMatx)

n_colors = length(break_seq)-1
break_seq_0 <- break_seq #save the originally specified break_seq to use in the color key creation
if (min(DifExpMatx)<break_seq[1]) {break_seq[1]=min(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
if (max(DifExpMatx)>break_seq[length(break_seq)]) {break_seq[length(break_seq)]=max(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
heat_map_colors <- colorRampPalette(c('blue','white','red'))(n_colors)
HeatmapName <- paste('heatmap',DistanceMethod,paste(ClusterMethod,'.pdf',sep=''),sep='_')
PDF_file <- paste(HeatmapDirectory,HeatmapName,sep='')
PdfW = 7
PdfH = 7
pdf(PDF_file,height=PdfH,width=PdfW) #not sure of the units of width and height
heatmap <- assemble_heatmap(GroupColorMatrix,DifExpMatx,colv,rowv,heat_map_colors,break_seq,label_rows,label_cols)
dev.off()
MakeColorKey(break_seq_0,heat_map_colors,HeatmapDirectory)
MakeGroupLegend(groups_corresponding,t(GroupColorMatrix),ColGroupsScheme,HeatmapDirectory)

return(list(C_col,C_row,groups_corresponding,DATA))
}
################################################################################################

#Function to return TRUE if the row does not have any NAs
NoNA <- function(vector)
{
  NoNA <- !is.na(sum(vector))
}

#Function to read a tab delimited text file into a data frame
read_txt_to_df <- function(txt_directory)
{
  DF <- read.table(file=txt_directory,head=TRUE,check.names=FALSE,sep='\t') #check.names=FALSE prevents an 'X' from being added to the numeric column names
  #Name the rows of the data frame as the genes given in the first column of the data frame
  RowNames <- as.character(DF[,1])
  rownames(DF) <- RowNames
  DF[,1] <- NULL #remove the first column of the data frame as it is no longer needed
  return(DF)
}

#Function to return row label_rows
FindRowLabels <- function(label_rows,DifExpMatx)
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

#Function to return label_cols
FindColLabels <- function(label_cols,DifExpMatx)
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
