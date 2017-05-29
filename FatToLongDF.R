FatToLongDF <- function(DATA,groups_corresponding)
{
    n_row_DATA <- nrow(DATA)
    n_col_DATA <- ncol(DATA)
    DATA_long_len <- nrow(DATA)*ncol(DATA)
    DATA_long <- data.frame(matrix(nrow=DATA_long_len,ncol=4))
    colnames(DATA_long) <- c('gene','sample','value','group')
    gene_names <- rownames(DATA)
    sample_names <- colnames(DATA)

    for (i in 1:n_col_DATA)
    {
        current_sample <- sample_names[i]
        start_iteration <- (i-1)*n_row_DATA + 1
        end_iteration <- i*n_row_DATA
        DATA_long[start_iteration:end_iteration,'gene'] <- rownames(DATA)
        DATA_long[start_iteration:end_iteration,'sample'] <- rep(current_sample,n_row_DATA)
        DATA_long[start_iteration:end_iteration,'value'] <- DATA[,current_sample]
        DATA_long[start_iteration:end_iteration,'group'] <- rep(groups_corresponding[current_sample,1],n_row_DATA)
    }

    return(DATA_long)
}
