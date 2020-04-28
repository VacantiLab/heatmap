RemoveMinValRows <- function(DATA,proportion)
{
    #This function removes rows where too many values are below 0 or too many values are above zero
    #The parapter proportion sets the limit of what fraction of numbers can be below or above zero
    #This is meant for log2 transformed ratios to a row mean or row median
    #proportion could be set to 0.7 or so

    n_samples <- length(colnames(DATA))
    rows_to_keep <- c()
    rows_discarded <- 0
    i = 0
    for (row in rownames(DATA))
    {
        if ((i >= 100)&(i %% 200 == 0))
        {
            print(paste('filtering ',toString(i),'th row',sep=''))
        }
        i = i+1
        n_less <- sum(DATA[row,] < 0)
        n_greater <- sum(DATA[row,] > 0)
        criteria <- (n_less/n_samples < proportion) & (n_greater/n_samples < proportion)
        if (criteria){rows_to_keep <- c(rows_to_keep,row)}
        if (!criteria){rows_discarded <- rows_discarded + 1}
    }
    DATA <- DATA[rows_to_keep,]
    print(paste(toString(rows_discarded),' rows discarded',sep=''))
    return(DATA)
}
