ArrangeData <- function(ColGroupsScheme,replicate_scheme,transformation,data,data_location,select_rows,select_groups,visualization,ddt,med_norm,handle_blanks)
# This function serves as a central data organization function for MakeVolcanoPlot, MakeBoxPlot, and MakeHeatMap
# med_norm specifies to median normalize columns
{
    #if more than one ColGroupsScheme is specified for a volcano plot because of a data dependent transformation (ddt)
    #    that needs to be disguised for the CheckStop
    ColGroupsScheme_holder = ColGroupsScheme
    if ((visualization == 'volcanoplot' |  visualization == 'boxplot') && !is.null(ddt) && length(ColGroupsScheme)>1)
    {
        ColGroupsScheme = ColGroupsScheme[ColGroupsScheme != ddt]
    }

    #Stop the program if the replicate scheme is in the ColGroupsScheme
    #Stop the program if more than one ColGroupsScheme or replicate_scheme is specified
    CheckStop(1,parameters=list(ColGroupsScheme,replicate_scheme,visualization))

    #Reset the ColGroupsScheme if it was disguised above
    ColGroupsScheme = ColGroupsScheme_holder

    #Include pertinent libraries
    library(ggplot2) #from ggplot2 package, allows the boxplot to be made

    #create the directory 'output' one level above the current directory for storing the heatmap
    output_directory <- StoreHeatmap()

    #specify if transformation occurs before or after groups of samples are excluded
    #    TRUE for after - default, FALSE for before
    transform_after_exclusion = TRUE

    #Input data
    if (!(class(data)=='data.frame')){data <- paste(data_location,'quantities.txt',sep='')}
    group_designations_file <- paste(data_location,'group_key.txt',sep='')
    group_color_designations_file <- paste(data_location,'group_color_key.txt',sep='')

    #Import the data and only keep selected rows and median-nomralize columns if specified
    print('opening data file')
    OpenDataFile_return <- OpenDataFile(data,select_rows,med_norm,handle_blanks)
    DATA <- OpenDataFile_return[[1]]
    select_rows <- OpenDataFile_return[[2]]

    #Retrieve the corresponding column groupings and keep only those specified
    #Also take the medians if there is a replicate scheme provided (this will not work if the replicate scheme is the only grouping scheme used)
    #There must be a ColGroupsScheme specified, as of now it cannot be FALSE

    #If the groups are packaged into lists, unpack them into a single array
        #This function does not affect the input if it is an array
    UnpackGroups_return <- UnpackGroups(select_groups)
    select_groups <- UnpackGroups_return[[1]] #this is now an array of the groups if it were input as a list
    group_divisions <- UnpackGroups_return[[2]] #this is the original input if it was a list, NULL if the input was not a list

    #Match the group names to the samples by referencing the group_key
    print('retrieving group and color assignments')
    ColGroupsScheme_concat <- c(ColGroupsScheme,replicate_scheme)
    RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme_concat,group_designations_file,group_color_designations_file,select_groups)
    groups_corresponding <- RetrieveGroups_return[[1]]
    GroupColorMatrix <- RetrieveGroups_return[[2]]
    COLOR_KEY <- RetrieveGroups_return[[3]]
    CheckStop(5,parameters=list(COLOR_KEY)) #makes sure each group name only has one color assignment

    #See if all of the specified input groups are actually specified in group_key.txt file
    CheckStop(2,parameters=list(select_groups,groups_corresponding))

    #transform if specified to do so before excluding groups or samples
    DATA_transformed_full <- NULL
    if (transform_after_exclusion == FALSE)
    {
        print('transforming data if necessary')
        DATA <- transform_data(DATA,transformation)
        DATA_transformed_full <- DATA
        gene_name <- rownames(DATA)
        n_gene <- length(gene_name)
    }


    #Select the groups that are considered for this box plot
    #    This is automatically performed on the first ColGroupsScheme provided
    print('selecting groups if necessary')
    SelectGroups_return <- SelectGroups(select_groups,DATA,ColGroupsScheme_concat,groups_corresponding,GroupColorMatrix,inclusion_grouping_scheme=ColGroupsScheme[1])
    #inclusion_grouping_scheme will need to be specified when more than one grouping scheme can be used such as in a heatmap
    DATA <- SelectGroups_return[[1]]
    groups_corresponding <- SelectGroups_return[[2]]
    GroupColorMatrix <- SelectGroups_return[[3]]


    #If you are concatonating groups, name the new groups and replace all of the groups they map to those with names
    #    There can only be one ColGroupsScheme provided if there is group concatonation
    #Also get corresponding colors for those new groups by taking the color that maps to the first sub-group of each concatonated group
    #This function does not affect the input if group_divisions is NULL (i.e. select_groups was not passed as a list to the original function)
    print('combining groups if necessary')
    ConcatonateGroups_return <- ConcatonateGroups(group_divisions,groups_corresponding,GroupColorMatrix,COLOR_KEY,concat_group_scheme=ColGroupsScheme)
    #concat_group_scheme has to be ColGroupsScheme because a volcano plot only has one ColGroupsScheme
    groups_corresponding <- ConcatonateGroups_return[[1]]
    GroupColorMatrix <- ConcatonateGroups_return[[2]]
    groups_concatonated <- ConcatonateGroups_return[[3]]
    colors_concatonated <- ConcatonateGroups_return[[4]]
    group_concationation <- is.list(group_divisions)


    #return median of replicates if specified to do so
    print('aggregating group members if necessary')
    MedianGroup_return <- MedianGroup(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme,ColGroupsScheme)
    DATA <- MedianGroup_return[[1]]
    groups_corresponding <- MedianGroup_return[[2]]
    GroupColorMatrix <- MedianGroup_return[[3]]


    #transform the columns based on other columns if specified to do so
    if (class(ddt)=='list')
    {
        print('transforming columns')
        TransformColumns_return <- TransformColumns(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme,ColGroupsScheme,ddt)
        DATA <- TransformColumns_return[[1]]
        groups_corresponding <- TransformColumns_return[[2]]
        GroupColorMatrix <- TransformColumns_return[[3]]
    }


    #median normalize within each of the groups as designated by the grouping scheme specified by the input ddt
    if (!is.null(ddt))
    {
        if (class(ddt) == 'character')
        {
            med_norm_scheme = ddt
            DATA = med_norm_within_groups(DATA,groups_corresponding,med_norm_scheme)
        }
    }

    #reset the ColGroupsScheme and groups_corresponding to what they would have been without a ddt
    #    may not be necessary
    if (visualization == 'volcanoplot' && !is.null(ddt) && length(ColGroupsScheme)>1)
    {
        ColGroupsScheme = ColGroupsScheme[ColGroupsScheme != ddt]
        cols_to_keep = colnames(groups_corresponding) != ddt
        groups_corresponding = groups_corresponding[,cols_to_keep,drop=FALSE]
    }

    variance_threshhold = NULL
    if (!is.null(variance_threshhold))
    {
        gene_vars <- apply(DATA,1,var)
        gene_means <- apply(DATA,1,mean)
        rel_var <- gene_vars/gene_means
        thresh_array <- rel_var > variance_threshhold
        DATA_rows_keep <- rownames(DATA)[thresh_array]
        DATA <- DATA[DATA_rows_keep,]
    }

    #transform if specified to do so after excluding groups or samples
    if (transform_after_exclusion == TRUE)
    {
        print('transforming data if necessary')
        DATA <- transform_data(DATA,transformation)
        gene_name <- rownames(DATA)
        n_gene <- length(gene_name)
    }

    #The rows with zero variance are removed because likely this causes NAs to be generated when calculating correlations
        #Most of the time this won't be an issue because all samples are not likely to have identical values for a measurement
        #However in measruing 10k protein abundances across few samples, this could happen by chance
    #print('searching for and removing rows with zero variance')
    #DATA <- remove_zero_var_rows(DATA)

    DATA_long <- NULL
    if (visualization=='boxplot' | visualization=='volcanoplot' | visualization=='scatterbar')
    {
        #put the data in long data frame format for plotting
        print('arranging data in long data frame format')
        DATA_long <- FatToLongDF(DATA,groups_corresponding)
    }

    FillColors = NULL #needs to be designated because it is returned
    group_order = NULL #need to be designated because it is returned

    if (visualization=='boxplot' | visualization=='volcanoplot')
    {
        #specify the order in which the groups will be plotted and ensure they map to their corresponding colors
        if (is.character(select_rows)){gene_name <- select_rows}
        OrderGroups_return <- OrderGroups(select_groups,group_concationation,groups_corresponding,GroupColorMatrix,COLOR_KEY,groups_concatonated,colors_concatonated,gene_name,DATA_long)
        DATA_long <- OrderGroups_return[[1]]
        FillColors <- OrderGroups_return[[2]]
        group_order <- OrderGroups_return[[3]]
    }

    sig_test_list <- list() #need to be designated because it is returned
    #only performed now if rows are not selected because if a gene is listed in select_rows that is not in the database, GetPs currently causes an error (2017-07-12)
    if (visualization=='boxplot' | visualization=='volcanoplot' | is.null(select_rows))
    {
        #calculate p-values: currently can only do if there are 2 groups with an equal-variance t-test
        #returns a data frame with the row naming the group pairwise comparison (e.g. basal-her2) and the column the parameter measured (e.g. glucose)
        print('performing significance tests')
        n_groups <- length(group_order)
        if (n_groups==2)
        {
            sig_test_list <- GetPs(group_order,gene_name,DATA,groups_corresponding,visualization)
        }
    }

    ArrangeData_return <- list(sig_test_list,output_directory,group_order,gene_name,DATA_long,FillColors,DATA,GroupColorMatrix,groups_corresponding,DATA_transformed_full)
    if (visualization == 'volcanoplot')
    {
        ArrangeData_return <- list(sig_test_list,output_directory,group_order,gene_name,DATA_long,FillColors,DATA,GroupColorMatrix,groups_corresponding,DATA_transformed_full,ColGroupsScheme)
    }
    return(ArrangeData_return)
}
