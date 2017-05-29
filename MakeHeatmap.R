MakeHeatMap <- function(data_location,ColGroupsScheme=FALSE,transformation=NULL,break_seq=NULL,replicate_scheme=NULL,DistanceMethod='pearson',ClusterMethod='ward.D2',data=NULL,select_rows=NULL,select_groups=NULL,label_rows=TRUE,label_cols=TRUE,rev_c_dend=FALSE)
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
#select_groups is a list of group designations that you want to consider in the heatmap (corresponding to the column names)
  #if NULL, all columns are considered
  #This is not operable yet
#rev_c_dend reverses the column dendrogram if TRUE
#label_rows specifies whether the rows should be labeled. If TRUE, all rows are labeled, if FALSE, no rows are labeled.
  #If it is a vector of strings, only rows with names within the vector of strings are labeled.
#This function also prints the heatmap as a pdf into the specified folder
#This function returns a list (dendrograms - or any name you specify) with the row and column dendrograms.
#The row dendrogram is called by dendrograms[[1]]

{
    #Include pertinent libraries
    library(gplots)
    library(heatmap.plus)

    #create the directory 'output' one level above the current directory for storing the heatmap
    HeatmapDirectory <- StoreHeatmap()

    #Input data
    if (!(class(data)=='data.frame')){data <- paste(data_location,'quantities.txt',sep='')}
    group_designations_file <- paste(data_location,'group_key.txt',sep='')
    group_color_designations_file <- paste(data_location,'group_color_key.txt',sep='')

    #Import the data and only keep selected rows if specified
    DATA <- OpenDataFile(data,select_rows)

    #Retrieve the corresponding column groupings and keep only those specified
    #Also take the medians if there is a replicate scheme provided (this may not work if the replicate scheme is the only grouping scheme used)
    RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme,group_designations_file,group_color_designations_file,select_groups,replicate_scheme)
    DATA <- RetrieveGroups_return[[1]]
    groups_corresponding <- RetrieveGroups_return[[2]]
    GroupColorMatrix <- RetrieveGroups_return[[3]]

    #Transform the data as specified
    DATA <- transform_data(DATA,transformation)

    #The rows with zero variance must be removed because likely because this causes NAs to be generated when calculating correlations
    #Most of the time this won't be an issue because all samples are not likely to have identical values for a measurement
        #However in measruing 10k protein abundances across few samples, this could happen by chance
    DATA <- remove_zero_var_rows(DATA)

    #Retrieve the clusters for creating the dendrograms
    ClusterData_return <- ClusterData(DATA,DistanceMethod,ClusterMethod,rev_c_dend)
    colv <- ClusterData_return[[1]]
    rowv <- ClusterData_return[[2]]
    C_col <- ClusterData_return[[3]]
    C_row <- ClusterData_return[[4]]

    #make numeric matrix out of data frame
    DifExpMatx=matrix(as.numeric(unlist(DATA)),nrow=nrow(DATA))
    rownames(DifExpMatx)<-rownames(DATA) #name the rows of the numeric matrix so they are transfered to the dendrogram in the heat map
    colnames(DifExpMatx)<-colnames(DATA) #name the columns of the numeric matrix so they are transfered to the dendrogram in the heat map

    #Get the row names and column names for the heatmap
    label_rows <- FindRowLabels(label_rows,DifExpMatx)
    label_cols <- FindColLabels(label_cols,DifExpMatx)

    #Draw and save the heatmap
    assemble_heatmap_return <- assemble_heatmap(GroupColorMatrix,DifExpMatx,colv,rowv,heat_map_colors,break_seq,label_rows,label_cols,HeatmapDirectory,DistanceMethod,ClusterMethod)
    heat_map_colors <- assemble_heatmap_return[[1]]
    dev.off() #turn off printing to the specified pdf

    #Draw and save the legends
    MakeColorKey(break_seq,heat_map_colors,HeatmapDirectory)
    MakeGroupLegend(groups_corresponding,t(GroupColorMatrix),ColGroupsScheme,HeatmapDirectory)

    #Return what could be used
    return(list(C_col,C_row,groups_corresponding,DATA))
}
