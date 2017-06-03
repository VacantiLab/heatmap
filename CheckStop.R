CheckStop <- function(check_flag)
{
    if (check_flag==1)
    {
        #Stop the program if the replicate scheme is in the ColGroupsScheme
        if (!is.null(ColGroupsScheme) && !is.null(replicate_scheme))
        {
            replicate_scheme_unique = sum(replicate_scheme %in% ColGroupsScheme)==0
            if(!replicate_scheme_unique){stop('custom message: The replicate_scheme is in the ColGroupsScheme. These should be mutually exclusive.')}
        }

        #Stop the program if more than one ColGroupsScheme or replicate_scheme is specified
        if (length(ColGroupsScheme)>1 || length(replicate_scheme)>1){stop('custom message: Cannot have more than one ColGroupsScheme or replicate_scheme for a boxplot')}
    }
}
