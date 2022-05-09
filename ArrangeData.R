ArrangeData <- function(ColGroupsScheme,
                       replicate_scheme,
                       transformation,
                       data,
                       data_location,
                       select_rows,
                       select_groups,
                       visualization,
                       ddt,
                       handle_blanks,
                       inclusion_grouping_scheme,
                       ttest,
                       select_rows_after_transform,
                       transform_after_column_exclusion,
                       FilterRowsForMeanSpread=FALSE,
                       GeneToCorrelate=NULL,
                       ratio_scheme_groups=NULL)

 # Inputs:
 # ColGroupsScheme (required): the name of the grouping scheme used indicates which column to take from group_key.txt
 #    There can only be one for the MakeBoxPlot() function
 # replicate_scheme (can be NULL): This specifies the grouping scheme that is used to specify groups of replicates
 #    This must NOT be a member of ColGroupsScheme, though it must be a grouping scheme defined in group_key.txt
 #        As such each member of this grouping scheme must also have colors specified in group_color_key.txt
 #    If this is specified, all members of a single group are treated as a single sample and the median values are used
 # transformation (can be NULL): This specifies how the rows should be transformed.
 #    Options include: 'log2', 'median_center_iqr_norm', and 'median_norm_log2_transform', and others.
 #    See transform_data.R for a complete list of possible transformations.
 # data (can be NULL): Is a data frame containing the data to be plotted if the data is passed as a data frame
 #       If this is used then the data is not passed through the quantities.txt file in data_location
 # data_location (required): a pathway to where the text file containing the data is stored, must have '/' at the end
 #    The data file must be named quantities.txt with the genes down the rows and sample names across the columns
 #    There must also be a group_key.txt file with the sample names down the rows and the grouping schemes across the columns
 #        The group names within each grouping scheme must not match any of those of another grouping scheme
 #        Each column must have a header
 #    There must also be a group_color_key.txt file with the group names down the rows and the colors beside them
 #        All of the group names are listed down the rows with no indication of their scheme membership (thus the names must be unique)
 #        Each column must have a header
 #    These files must all be tab delimited
 # select_rows (can be NULL): These are the rows that you want to select before transforming - the transformation will be done on these rows
 #    If it is NULL, all genes included
 #    If it is a path to a .txt file, that file must contain the gene names, one on each line, and they will be used as the select_rows
 #    The transformations will not be done on the whole data set and this does impact the values
 # select_groups (can be NULL): This can be an array of the group names that are to be plotted or a list of arrays of group names or NULL
 #    If it is an array of group names, those groups are the only ones plotted
 #    If it is a list of arrays of group names, groups in the same array are combined into a single group
 #    If it is NULL, all groups in the ColGroupsScheme are plotted
 #    These can be groups outside of ColGroupsScheme, but the scheme must then be specified as the inclusion_grouping_scheme
 #    The spcification of transform_after_column_exclusion determines if this impacts data transformation/normalization
 # visualization (required): A string describing the type of plot
 #    Can be one of the visualizations this function is used to produce such as: 'volcanoplot'
 # ddt (can be NULL): data dependent transformation; this is a grouping scheme that all samples within that group are normalized to its median and then log2 transformed
 #    It must be one of the groups specified in the ColGroupsScheme
 #    It is not presented as a ColGroupsScheme, it is just used for normalization purposes
 #        For example, say you have cell lines control and treated
 #            You can specify to normalize within cell lines and then use select_groups and inclusion_grouping_scheme to plot only the treated samples
 #                The result would be the treatment response for each cell line
 #    transform_after_column_exclusion must be FALSE because the transformation would also occur after DDT which doesn't make sense
 #        Sample loading should be accounted for before DDT
 #        The transformation should also be linear because DDT log2 transforms resulting ratios
 # handle_blanks (required): Used in OpenDataFile in the SupportingFunctions.R file.
 #                It is a string that specifies what to do when blank values are encountered
 #                Can be 'remove_row' or 'replace_with_rowmin'
 # inclusion_grouping_scheme (can be NULL): This is the grouping scheme that you want to specify to select the columns with
 #    It can remain as the default NULL and the groups will be selected based on the first grouping scheme in ColGroupsScheme if select_groups is specified
 #    It can be outside of the ColGroupsScheme as well
 # ttest (required): Indicates whether ttests should be performed by being either TRUE of FALSE
 #    only works if there are two select_groups
 # select_rows_after_transform (can be NULL): these rows are selected to be plotted after tranformations have been completed
 #     This does not impact values
 # transform_after_column_exclusion (required): if TRUE, transformations happen after groups (columns) are selected, not on the whole data set
 #     If TRUE, select_groups impacts the values
 #     If FALSE, select_groups does not impact the values
 #     If TRUE, this may impact values if ddt is specified
 #     If FALSE, this will not impact values if ddt is specified
 # FilterRowsForMeanSpread (can be NULL): If it is a decimal fraction (such as 0.7), rows are filtered out who have more than 70% of their values below or above 0
 #     This is meant for log2 transformed data where there may be a small cluster of columns driving the profile
 #     Corrects red or blue streaks across the heat map

 #Outputs
 # The Data frame retunred by this function is the 7th output of the returned list
 #    Thus in the calll: ADR <- ArrangeData()
 #        DATA <- ADR[[7]]


# This function serves as a central data organization function for MakeVolcanoPlot, MakeBoxPlot, and MakeHeatMap
{
    #if more than one ColGroupsScheme is specified for a volcano plot because of a data dependent transformation (ddt)
    #    that needs to be disguised for the CheckStop
    ColGroupsScheme_holder = ColGroupsScheme
    if ((visualization == 'volcanoplot' |  visualization == 'boxplot' | visualization == 'ScatterLinePlot') && !is.null(ddt) && length(ColGroupsScheme)>1)
    {
        ColGroupsScheme = ColGroupsScheme[ColGroupsScheme != ddt]
    }

    #Stop the program if the replicate scheme is in the ColGroupsScheme
    #Stop the program if more than one ColGroupsScheme or replicate_scheme is specified
    CheckStop(1,parameters=list(ColGroupsScheme,replicate_scheme,visualization,inclusion_grouping_scheme))

    #Reset the ColGroupsScheme if it was disguised above
    ColGroupsScheme = ColGroupsScheme_holder

    #Include pertinent libraries
    library(ggplot2) #from ggplot2 package, allows the boxplot to be made

    #create the directory 'output' one level above the current directory for storing the heatmap
    output_directory <- StoreHeatmap()

    #specify if transformation occurs before or after groups of samples are excluded
    #    TRUE for after - default, FALSE for before
    #transform_after_column_exclusion = TRUE

    #Input data
    if (!(class(data)=='data.frame')){data <- paste(data_location,'quantities.txt',sep='')}
    group_designations_file <- paste(data_location,'group_key.txt',sep='')
    group_color_designations_file <- paste(data_location,'group_color_key.txt',sep='')

    #Import the data and only keep selected rows and median-nomralize columns if specified
    print('opening data file')
    OpenDataFile_return <- OpenDataFile(data,select_rows,handle_blanks)
    DATA <- OpenDataFile_return[[1]]
    select_rows <- OpenDataFile_return[[2]]
    DATA_original <- OpenDataFile_return[[3]]

    #Add the data you want to co-analyze if doing an edge_table for a network analysis
    #This is hardcoded for now
    # if(visualization=='edge_table')
    # {
    #     DATA['m3malate',c('bt20_rotenone','hcc1419_rotenone','mcf7_rotenone','t47d_rotenone','mdamb157_rotenone','bt20_control','hcc1419_control','mcf7_control','t47d_control','mdamb157_control')] <- c(0.11458,0.05459,0.09237,0.019732,0.21588,0.061236,0.094963,0.179303,0.065843,0.171299)
    # }

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

    #The replicate scheme needs to be considered to find the medians
    #    It is not considered for plotting different groups
    ColGroupsScheme_concat <- c(ColGroupsScheme,replicate_scheme)

    #The inclusion grouping scheme needs to be added to ColGroupsScheme if it is not already there
    #    It's removed later
    if (!(inclusion_grouping_scheme %in% ColGroupsScheme) && !is.null(inclusion_grouping_scheme))
    {
        ColGroupsScheme_concat <- c(ColGroupsScheme,inclusion_grouping_scheme)
    }

    #Column grouping schemes needed to consider will include the inclusion_grouping_scheme as well
    #    However this only needs to be included if it is different from the replicate scheme (otherwise it is already there from above)
    if (!(is.null(replicate_scheme)) & !(is.null(inclusion_grouping_scheme))) #logical statements return emtpy values if they deal with NULL
    {
        if (replicate_scheme != inclusion_grouping_scheme){ColGroupsScheme_concat <- c(ColGroupsScheme,inclusion_grouping_scheme)}
    }

    RetrieveGroups_return <- RetrieveGroups(DATA,ColGroupsScheme_concat,group_designations_file,group_color_designations_file,select_groups)
    groups_corresponding <- RetrieveGroups_return[[1]]
    GroupColorMatrix <- RetrieveGroups_return[[2]]
    COLOR_KEY <- RetrieveGroups_return[[3]]
    CheckStop(5,parameters=list(COLOR_KEY)) #makes sure each group name only has one color assignment

    #See if all of the specified input groups are actually specified in group_key.txt file
    CheckStop(2,parameters=list(select_groups,groups_corresponding,inclusion_grouping_scheme))
    
    #Filter out rows depending on their distribution around the mean
    DATA_after_filter <- NULL
    if (class(FilterRowsForMeanSpread) == 'numeric')
    {
        proportion <- FilterRowsForMeanSpread
        DATA <- RemoveMinValRows(DATA,proportion)
        DATA_after_filter <- DATA
    }

    #transform if specified to do so before excluding groups or samples
    DATA_transformed_full <- NULL
    if (handle_blanks == 'remove_row_after_group_exclusion'){transform_after_column_exclusion == TRUE}
    #    You cannot transform data if there are still NAs in the data frame
    if (transform_after_column_exclusion == FALSE)
    {
        print('transforming data if necessary')
        DATA <- transform_data(DATA,transformation,select_rows_after_transform)
        DATA_transformed_full <- DATA
        gene_name <- rownames(DATA)
        n_gene <- length(gene_name)
    }

    #Perform the data dependent transformation if specified to do so
    #This is done before any columns are excluded
    #    There may never be a need to do it after columnms are excluded because the values only depend on other values within the group
    #    This should be done before the exclusion and after the transformation. The transformation should be linear because DDT log2 transforms the ratios returned

    #Make sure the trasnformation occurs before column exclusion (? 20210407 moved ddt to after column exclusion) and the transformation is linear because DDT includes log2 transformation of returned ratios
    CheckStop(8,parameters=list(ddt,transformation,transform_after_column_exclusion))
  

    #Select the groups that are considered for this box plot
    #    This is automatically performed on the first ColGroupsScheme provided
    if (is.null(inclusion_grouping_scheme)){inclusion_grouping_scheme=ColGroupsScheme[1]}
    print('selecting groups if necessary')
    SelectGroups_return <- SelectGroups(select_groups,DATA,ColGroupsScheme_concat,groups_corresponding,GroupColorMatrix,inclusion_grouping_scheme,visualization)
    #inclusion_grouping_scheme will need to be specified when more than one grouping scheme can be used such as in a heatmap
    DATA <- SelectGroups_return[[1]]
    groups_corresponding <- SelectGroups_return[[2]]

    
    # Remove rows with NA entries after the specified groups have been excluded
    #     If it is specified to do so
    #     Note in this case the transformation must occur after groups have been excluded
    #         This will be specified automatically
    if (handle_blanks == 'remove_row_after_group_exclusion')
    {
      has_no_na_row_indices <- apply(DATA,1,NoNA)
      DATA <- DATA[has_no_na_row_indices,,drop=FALSE]
    }
    
    if (!is.null(ddt))
    {
      PerformDDT_return <- PerformDDT(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme,ColGroupsScheme,ddt,ratio_scheme_groups=ratio_scheme_groups)
      DATA <- PerformDDT_return[[1]]
      groups_corresponding <- PerformDDT_return[[2]]
      GroupColorMatrix <- PerformDDT_return[[3]]
    }
    

    #remove the inclusion grouping scheme as a category in groups_corresponding if it is not part of the ColGroupsScheme
    #    this means your plot does not consider the inclusion grouping scheme, it was just used to select the data to plot
    #    i.e. I only want to plot data of a certain cell type, but I do not want cell type to be a category in my plot
    GroupColorMatrix <- SelectGroups_return[[3]]
    if (!(inclusion_grouping_scheme %in% ColGroupsScheme) && ((replicate_scheme != inclusion_grouping_scheme) || is.null(replicate_scheme)))
    {
      current_group_categories <- colnames(groups_corresponding)
      group_categories_to_keep <- current_group_categories != inclusion_grouping_scheme
      # do not consider the inclusion grouping scheme as a group to plot
      groups_corresponding <- groups_corresponding[,group_categories_to_keep,drop=FALSE]
      GroupColorMatrix <- GroupColorMatrix[,group_categories_to_keep,drop=FALSE]
    }

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
    if (transform_after_column_exclusion == TRUE)
    {
        print('transforming data if necessary')
        DATA <- transform_data(DATA,transformation,select_rows_after_transform)
        gene_name <- rownames(DATA)
        n_gene <- length(gene_name)
    }

    #The rows with zero variance are removed because likely this causes NAs to be generated when calculating correlations
        #Most of the time this won't be an issue because all samples are not likely to have identical values for a measurement
        #However in measruing 10k protein abundances across few samples, this could happen by chance
    #print('searching for and removing rows with zero variance')
    #DATA <- remove_zero_var_rows(DATA)


    DATA_long <- NULL
    if (visualization=='boxplot' | visualization=='volcanoplot' | visualization=='scatterbar' | visualization=='ScatterLinePlot')
    {
        #put the data in long data frame format for plotting
        print('arranging data in long data frame format')
        DATA_long <- FatToLongDF(DATA,groups_corresponding,ddt)
    }

    FillColors = NULL #needs to be designated because it is returned
    group_order = NULL #need to be designated because it is returned

    if (visualization=='boxplot' | visualization=='volcanoplot')
    {
        #specify the order in which the groups will be plotted and ensure they map to their corresponding colors
        #if (is.character(select_rows)){gene_name <- select_rows}
        groups_to_order <- select_groups
        #if you are not selecting groups based on the ColGroupsScheme, groups_to_order should be the default NULL (it will be assumed to be the order of the ColGroupsScheme?)
        if (!(inclusion_grouping_scheme %in% ColGroupsScheme)){groups_to_order <- NULL}
        groups_corresponding_ForOrder <- groups_corresponding
        #if there is a ddt, it needs to be removed from the groups_corresponding when determining the order of the groups because ddt is not actually a group but a reference for a transformation
        if (!is.null(ddt) & visualization=='heatmap'){groups_corresponding_ForOrder <- groups_corresponding[,colnames(groups_corresponding) != ddt,drop=FALSE]}
        OrderGroups_return <- OrderGroups(groups_to_order,group_concationation,groups_corresponding_ForOrder,GroupColorMatrix,COLOR_KEY,groups_concatonated,colors_concatonated,gene_name,DATA_long,ddt,visualization)
        DATA_long <- OrderGroups_return[[1]]
        FillColors <- OrderGroups_return[[2]]
        group_order <- OrderGroups_return[[3]]
    }


    sig_test_list <- list() #need to be designated because it is returned

    #only shoule be performed now if rows are not selected because if a gene is listed in select_rows that is not in the database, GetPs currently causes an error (2017-07-12)
    if ((visualization=='boxplot' | visualization=='volcanoplot' | is.null(select_rows)) & ttest)
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
    
    #only shoule be performed now if rows are not selected because if a gene is listed in select_rows that is not in the database, GetPs currently causes an error (2017-07-12)
    if (visualization=='CorrelationVolcanoPlot')
    {
      print(paste('calculating pairwise correlations with ',GeneToCorrelate))
      {
        Correlations <- PerformCorrelations(DATA,GeneToCorrelate)
      }
    }
    
    ArrangeData_return <- list(sig_test_list,
                               output_directory,
                               group_order,
                               gene_name,
                               DATA_long,
                               FillColors,
                               DATA,
                               GroupColorMatrix,
                               groups_corresponding,
                               DATA_transformed_full,
                               DATA_original)
    
    if (visualization == 'volcanoplot')
    {
        ArrangeData_return <- list(sig_test_list,
                                   output_directory,
                                   group_order,
                                   gene_name,
                                   DATA_long,
                                   FillColors,
                                   DATA,
                                   GroupColorMatrix,
                                   groups_corresponding,
                                   DATA_transformed_full,
                                   ColGroupsScheme,
                                   DATA_after_filter)
    }
    
    if (visualization == 'CorrelationVolcanoPlot')
    {
      ArrangeData_return <- list(Correlations,
                                 output_directory,
                                 DATA)
    }
    return(ArrangeData_return)
}

###############################################################################

NoNA <- function(vector)
  #Function to return TRUE if the row does not have any NAs
{
  NoNA <- !is.na(sum(vector))
}
