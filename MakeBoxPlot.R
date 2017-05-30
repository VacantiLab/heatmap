MakeBoxPlot <- function(data_location,ColGroupsScheme=FALSE,transformation,data=NULL,select_rows=NULL,select_groups=NULL,replicate_scheme=NULL)
# data_file: A text file where the rows are patients and the columns are genes, includes the directory of where the file is located
# box_plot_directory: The location of where the output box plot will be saved, has a "/" at the end
# group_designations_file: A tab delimited text file where each entry corresponds to the group each column in the data_file belongs to
# SelectGenes: A list of the genes, as strings, that are plotted along the x-axis. The default value of 'all' will plot every gene in data_file
# TakeLog2: If TRUE all expression values are log2 transformed, if false values are not log2 transformed
{
    #Include pertinent libraries
    library(ggplot2) #from ggplot2 package
    library(tidyr) #from tidyr package

    #create the directory 'output' one level above the current directory for storing the heatmap
    BoxDirectory <- StoreHeatmap()

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
    groups_corresponding <- t(RetrieveGroups_return[[2]])
    GroupColorMatrix <- RetrieveGroups_return[[3]]
    COLOR_KEY <- RetrieveGroups_return[[4]]

    #Transform the data as specified
    DATA <- transform_data(DATA,transformation)

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
    DATA_long$group <- factor(DATA_long$group,group_order_original) #this sets the order of the groups to match group_order_original

    #Make the plot
    b <- assemble_box_plot(DATA_long,FillColors,BoxDirectory)

    return(b)
}
