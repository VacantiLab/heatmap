MakeBoxPlot <- function(data_location,ColGroupsScheme=FALSE,transformation,data=NULL,select_rows=NULL,select_groups=NULL,replicate_scheme=NULL)
# data_location: a pathway to where the text file containing the data is stored
#    The data file must be named quantities.txt with the genes down the rows and sample names across the columns
#    There must also be a group_key.txt file with the sample names across the columns and the grouping schemes down the rowss
#        The group names within each grouping scheme must not match any of those of another grouping scheme
#    There must also be a group_color_key.txt file with the group names across the columns and the colors beneath them
#        All of the group names are listed across the column headers with no indication of their scheme membership (thus the names must be unique)
# box_plot_directory: The location of where the output box plot will be saved, has a "/" at the end
# ColGroupsScheme: the name of the grouping scheme used, indicates which row to take from group_key.txt
#    This needs to be specified, currently this cannot be run without a grouping scheme
# transformation: This specifies how the rows should be transformed.
#    Options include: 'log2', 'median_center_iqr_norm', and 'median_norm_log2_transform'
# data: Is a data frame containing the data to be plotted if the data is passed as a data frame instead of through the quantities.txt file in data_location
# select_rows: These are the genes that you want to plot
#    #If it is NULL, all genes are plotted
# select_groups: This can be an array of the group names (from the provided ColGroupsScheme) that are to be plotted or a list of arrays of group names or NULL
#    If it is an array of group names, those groups are the only ones plotted
#    If it is a list of arrays of group names, groups in the same array are combined into a single group
#    If it is NULL, all groups in the ColGroupsScheme are plotted
# replicate_scheme: This specifies the ColGroupsScheme that is used to specify groups of replicates
#    More than one A second ColGroupsScheme must be specified in addition to the replicate_scheme
#    If this is specified, all members of a single group are treated as a single sample and the median values are used
#    This is untested for MakeBoxPlot
{
    #Include pertinent libraries
    library(ggplot2) #from ggplot2 package, allows the boxplot to be made

    #create the directory 'output' one level above the current directory for storing the heatmap
    BoxDirectory <- StoreHeatmap()

    #Input data
    if (!(class(data)=='data.frame')){data <- paste(data_location,'quantities.txt',sep='')}
    group_designations_file <- paste(data_location,'group_key.txt',sep='')
    group_color_designations_file <- paste(data_location,'group_color_key.txt',sep='')

    #Import the data and only keep selected rows if specified
    DATA <- OpenDataFile(data,select_rows)

    #Retrieve the corresponding column groupings and keep only those specified
    #Also take the medians if there is a replicate scheme provided (this will not work if the replicate scheme is the only grouping scheme used)
    #There must be a ColGroupsScheme specified, as of now it cannot be FALSE
    group_concatonation = is.list(select_groups)
    group_divisions = NULL

    #Unpack the group names if they are packed in a list of arrays to be concatonated
    #If you are concatonating groups, you need to have an un-concatonated array of them to match to the group_key
    if (group_concatonation)
    {
        group_divisions = select_groups
        select_groups = do.call(c,select_groups)
    }

    #Match the group names to the samples by referencing the group_key
    RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme,group_designations_file,group_color_designations_file,select_groups,replicate_scheme)
    DATA <- RetrieveGroups_return[[1]]
    groups_corresponding <- t(RetrieveGroups_return[[2]])
    GroupColorMatrix <- RetrieveGroups_return[[3]]
    COLOR_KEY <- RetrieveGroups_return[[4]]

    #If you are concatonating groups, name the new groups and replace all of the groups they map to with those names
    #Also get corresponding colors for those new groups by taking the color that maps to the first sub-group of each concatonated group
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

    #Transform the data as specified
    DATA <- transform_data(DATA,transformation)
    gene_name <- rownames(DATA)
    n_gene <- length(gene_name)

    #put the data in long data frame format for plotting
    DATA_long <- FatToLongDF(DATA,groups_corresponding)

    #specify the order in which the groups will be plotted and ensure they map to their corresponding colors
    group_order <- matrix(as.character(unique(groups_corresponding)),ncol=1)
    FillColors <- matrix(as.character(COLOR_KEY[group_order]),ncol=1)

    #If groups to be plotted are specified, their order specifies the order they will be plotted in
    #This allows the user to control the order groups are plotted in
    if (is.character(select_groups))
    {
        group_order <- matrix(as.character(select_groups),ncol=1)
        FillColors <- matrix(as.character(COLOR_KEY[select_groups]),ncol=1)
    }

    #The group list and corresponding FillColors for ordering must reflect the combined groups if there are combined groups
    #If there are combined groups, their order of plotting cannot be specifed (yet) and it follows from where their members appear in the data
    if (is.list(group_divisions))
    {
        group_order <- unique(groups_corresponding)
        FillColors <- as.matrix(unique(GroupColorMatrix[,1]),ncol=1)
    }

    #Implement the group ordering
    DATA_long$group <- factor(DATA_long$group,group_order) #this sets the order of the groups to match group_order
    DATA_long$gene <- factor(DATA_long$gene,gene_name) #this sets the order of the genes to match gene_names

    #calculate p-values: currently can only do if there are 2 groups with an equal-variance t-test
    #returns a data frame with the row naming the group pairwise comparison (e.g. basal-her2) and the column the parameter measured (e.g. glucose)
    n_groups <- length(group_order)
    p_val_df <- data.frame(matrix(nrow=1,ncol=n_gene))
    if (n_groups==2)
    {
        p_val_df <- GetPs(group_order,n_gene,gene_name,DATA_long)
    }

    #Find the y-limits for the boxplot based on the data
    y_bounds <- get_y_bounds(group_order,gene_name,DATA_long)

    #Make the plot
    b <- assemble_box_plot(DATA_long,FillColors,BoxDirectory,y_bounds)

    #assemble variables to return
    MakeBoxPlot_return <- list(p_val_df)

    return(MakeBoxPlot_return)
}
