MakeBoxPlot <- function(data_location,ColGroupsScheme=FALSE,transformation,data=NULL,select_rows=NULL,select_groups=NULL,replicate_scheme=NULL)
# data_location: a pathway to where the text file containing the data is stored
#    The data file must be named quantities.txt with the genes down the rows and sample names across the columns
#    There must also be a group_key.txt file with the sample names across the columns and the grouping schemes down the rowss
#        The group names within each grouping scheme must not match any of those of another grouping scheme
#    There must also be a group_color_key.txt file with the group names across the columns and the colors beneath them
#        All of the group names are listed across the column headers with no indication of their scheme membership (thus the names must be unique)
# box_plot_directory: The location of where the output box plot will be saved, has a "/" at the end
# ColGroupsScheme: the name of the grouping scheme used, indicates which row to take from group_key.txt
# transformation:
# data:
# select_rows:
# select_groups: This can be an array of the group names (from the provided ColGroupsScheme) that are to be plotted or a list of arrays of group names or NULL
#    
# replicate_scheme:
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
    if (group_concatonation)
    {
        group_divisions = select_groups
        select_groups = do.call(c,select_groups)
    }
    RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme,group_designations_file,group_color_designations_file,select_groups,replicate_scheme)
    DATA <- RetrieveGroups_return[[1]]
    groups_corresponding <- t(RetrieveGroups_return[[2]])
    GroupColorMatrix <- RetrieveGroups_return[[3]]
    COLOR_KEY <- RetrieveGroups_return[[4]]

    groups_concatonated <- lapply(group_divisions,paste,collapse=':')
    groups_concatonated <- unlist(groups_concatonated)
    n_groups_concatonated <- length(groups_concatonated)
    for (i in 1:n_groups_concatonated)
    {
        concatonate_indices <- groups_corresponding %in% group_divisions[[i]]
        groups_corresponding[concatonate_indices] <- groups_concatonated[i]
        GroupColorMatrix[concatonate_indices,1] <- COLOR_KEY[1,group_divisions[[i]][1]]
    }

    #Transform the data as specified
    DATA <- transform_data(DATA,transformation)
    gene_name <- rownames(DATA)
    n_gene <- length(gene_name)

    #put the data in long data frame format for plotting
    DATA_long <- FatToLongDF(DATA,groups_corresponding)


    group_order_original_plus <- colnames(COLOR_KEY)
    groups <- unique(groups_corresponding)
    group_order_indices <- match(groups,group_order_original_plus)
    group_order_indices_sorted <- sort(group_order_indices)
    group_order_original <- group_order_original_plus[group_order_indices_sorted]
    if (is.character(select_groups))
    {
        indices_to_keep = group_order_original %in% select_groups
        group_order_original = group_order_original[indices_to_keep]
    }

    FillColors <- matrix(as.character(COLOR_KEY[group_order_original]),ncol=1)

    #the group list and corresponding FillColors for ordering must reflect the combined groups if there are combined groups
    if (is.list(group_divisions))
    {
        group_order_original <- unique(groups_corresponding)
        FillColors <- as.matrix(unique(GroupColorMatrix[,1]),ncol=1)
    }


    DATA_long$group <- factor(DATA_long$group,group_order_original) #this sets the order of the groups to match group_order_original
    DATA_long$gene <- factor(DATA_long$gene,gene_name)


    #calculate p-values: currently can only do if there are 2 groups
    #returns a data frame with the row naming the group pairwise comparison (e.g. basal-her2) and the column the parameter measured (e.g. glucose)
    n_groups <- length(group_order_original)
    p_val_df <- data.frame(matrix(nrow=1,ncol=n_gene))
    if (n_groups==2)
    {
        p_val_df <- GetPs(group_order_original,n_gene,gene_name,DATA_long)
    }

    #Find the y-limits for the boxplot based on the data
    y_bounds <- get_y_bounds(group_order_original,gene_name,DATA_long)

    #Make the plot
    b <- assemble_box_plot(DATA_long,FillColors,BoxDirectory,y_bounds)

    #assemble variables to return
    MakeBoxPlot_return <- list(p_val_df)

    return(MakeBoxPlot_return)
}
