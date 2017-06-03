CheckStop <- function(check_flag,parameters)
{
    if (check_flag==1)
    {
        #Stop the program if the replicate scheme is in the ColGroupsScheme
        ColGroupScheme = parameters[[1]]
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
    if (check_flag==2)
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
    if (check_flag==3)
    {
        GROUP_KEY = parameters[[1]]
        ColGroupsScheme = parameters[[2]]
        possible_group_schemes <- rownames(GROUP_KEY)
        all_ColGroupScheme_real <- CheckAllIn(ColGroupsScheme,possible_group_schemes)
        if (!all_ColGroupScheme_real){stop('Custom Message: a specified ColGroupsScheme or the replicate_scheme does not exist.')}
    }

    #Make sure all of the specified select_groups are actually in the grouping scheme specified to select the groups from
    if (check_flag==4)
    {
        select_groups = parameters[[1]]
        groups_corresponding = parameters[[2]]
        inclusion_grouping_scheme = parameters[[3]]
        select_groups_in_same_scheme <- CheckAllIn(select_groups,groups_corresponding[,inclusion_grouping_scheme])
        if (!select_groups_in_same_scheme){stop('Custom Message: Not all select_groups are in the same ColGroupsScheme. You cannot specify to consider only group members from multiple grouping schemes.')}
    }

    #Make sure each group only has one color assignment
    if (check_flag==5)
    {
        COLOR_KEY = parameters[[1]]
        group_names <- colnames(COLOR_KEY)
        unique_group_names <- unique(group_names)
        if (length(group_names)!=length(unique_group_names)){stop('Custom Message: One or more groups have more than one color assignment.')}
    }
}
