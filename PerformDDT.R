PerformDDT <- function(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme,ColGroupsScheme,ddt)
{

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

    PerformDDT_return <- list(DATA,groups_corresponding,GroupColorMatrix)
    return(PerformDDT_return)


}
