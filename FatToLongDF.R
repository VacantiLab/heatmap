FatToLongDF <- function(DATA,groups_corresponding,ddt)
{
    n_row_DATA <- nrow(DATA)
    n_col_DATA <- ncol(DATA)
    DATA_long_len <- nrow(DATA)*ncol(DATA)
    DATA_long <- data.frame(matrix(nrow=DATA_long_len,ncol=3))
    colnames(DATA_long) <- c('gene','sample','value')
    gene_names <- rownames(DATA)
    sample_names <- colnames(DATA)
    not_ddt_group <- colnames(groups_corresponding)[1]

    if (!is.null(ddt))
    {
      not_ddt_group <- colnames(groups_corresponding)[colnames(groups_corresponding)!=ddt]
    }

    for (i in 1:n_col_DATA)
    {
        if (i%%100 == 0){print(paste('converting ',toString(i),'th column to long data frame',sep=''))}
        current_sample <- sample_names[i]
        start_iteration <- (i-1)*n_row_DATA + 1
        end_iteration <- i*n_row_DATA
        DATA_long[start_iteration:end_iteration,'gene'] <- rownames(DATA)
        DATA_long[start_iteration:end_iteration,'sample'] <- rep(current_sample,n_row_DATA)
        DATA_long[start_iteration:end_iteration,'value'] <- DATA[,current_sample]
        if (!is.null(groups_corresponding))
        {
            DATA_long[start_iteration:end_iteration,'group'] <- rep(groups_corresponding[current_sample,not_ddt_group],n_row_DATA)
            if (!is.null(ddt)){DATA_long[start_iteration:end_iteration,'ddt_group'] <- rep(groups_corresponding[current_sample,ddt],n_row_DATA)}
        }
    }
    return(DATA_long)
}
