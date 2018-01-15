CheckStop <- function(check_flag,parameters)
#Stop the program and display an informative error message if there are inconsistencies with the input parameters
{
    if (check_flag==1) #In MakeBoxPlot()
    {
        #Stop the program if the replicate scheme is in the ColGroupsScheme
        ColGroupsScheme = parameters[[1]]
        replicate_scheme = parameters[[2]]
        if (!is.null(ColGroupsScheme) && !is.null(replicate_scheme))
        {
            replicate_scheme_unique = sum(replicate_scheme %in% ColGroupsScheme)==0
            if(!replicate_scheme_unique){stop('Custom Message: The replicate_scheme is in the ColGroupsScheme. These should be mutually exclusive.')}
        }

        #Stop the program if more than one ColGroupsScheme or replicate_scheme is specified
        if (length(ColGroupsScheme)>1 || length(replicate_scheme)>1){stop('Custom Message: Cannot have more than one ColGroupsScheme or replicate_scheme for a boxplot')}

    }

    #See if all of the specified input groups are actually specified in group_key.txt file
    if (check_flag==2) #In MakeBoxPlot
    {
        select_groups = parameters[[1]]
        groups_corresponding = parameters[[2]]
        if (!is.null(select_groups) && !is.null(groups_corresponding))
        {
            all_select_groups_used <- sum(select_groups %in% groups_corresponding) == length(select_groups)
            if (!all_select_groups_used){stop('Custom Message: select_groups contains a member that no sample is assinged to in the specified ColGroupsScheme - check spelling')}
        }
    }

    #Make sure the specified ColGroupsScheme entries and replicate_scheme entries all exist
    #Ensure all sample names in the quantites.txt file map to a grouyp in the group_key.txt file
    if (check_flag==3) #In RetrieveGroups()
    {
        GROUP_KEY = parameters[[1]]
        ColGroupsScheme = parameters[[2]] #this includes the replicate_scheme because it was concatonated previously
        DATA = parameters[[3]]
        possible_group_schemes <- colnames(GROUP_KEY)
        all_ColGroupsScheme_real <- CheckAllIn(ColGroupsScheme,possible_group_schemes)
        if (!all_ColGroupsScheme_real){stop('Custom Message: a specified ColGroupsScheme or the replicate_scheme does not exist.')}
        all_samples_map_to_group <- CheckAllIn(colnames(DATA),rownames(GROUP_KEY))
        if (!all_samples_map_to_group){stop('Custom Message: Not all of the samples map to a group, check sample naming consistency in the quantities.txt and group_key.txt files')}
    }

    #Make sure all of the specified select_groups are actually in the grouping scheme specified to select the groups from
    if (check_flag==4) #In SelectGroups()
    {
        select_groups = parameters[[1]]
        groups_corresponding = parameters[[2]]
        inclusion_grouping_scheme = parameters[[3]]
        select_groups_in_same_scheme <- CheckAllIn(select_groups,groups_corresponding[,inclusion_grouping_scheme])
        if (!select_groups_in_same_scheme){stop('Custom Message: Not all select_groups are in the same ColGroupsScheme. You cannot specify to consider only group members from multiple grouping schemes.')}
    }

    #Make sure each group only has one color assignment
    if (check_flag==5) #In MakeBoxPlot()
    {
        COLOR_KEY = parameters[[1]]
        group_names <- rownames(COLOR_KEY)
        unique_group_names <- unique(group_names)
        if (length(group_names)!=length(unique_group_names)){stop('Custom Message: One or more groups have more than one color assignment.')}
    }

    #Ensure the sample->group and group->color mappings have consistent group names
    #The replicate_scheme is considered a grouping scheme here because this is after it is appended
    #Thus each member the replicate_scheme must have a color assigned to it in the group_color_key.txt file
    if (check_flag==6) #In GetGroupColorList()
    {
        COLOR_KEY <- parameters[[1]]
        GROUP_KEY <- parameters[[2]]
        ColGroupsScheme <- parameters[[3]]
        group_names_from_color_key <- rownames(COLOR_KEY)
        group_names_from_group_key <- as.matrix(GROUP_KEY[,ColGroupsScheme],nrow=length(ColGroupsScheme))
        group_names_consistent <- CheckAllIn(group_names_from_group_key,group_names_from_color_key)
        if (!group_names_consistent)
        {
            print(COLOR_KEY)
            print(GROUP_KEY)
            stop('The group names in group_key.txt are not consistent with the group names in the group_color_key.txt.')
        }
    }

    if (check_flag==7) #In MakeVolcanoPlot
    {
        transformation <- parameters[[1]]
        if (is.character(transformation))
        {
            if (transformation == 'log2'){stop('Custom Message: The transformation cannot be log2 for a volcano plot.')}
        }
    }
}
