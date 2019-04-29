CommonRows <- function(df1,df2)
# note this is not needed for any plotting Function
#     it is written only for convenience to make mRNA and protein plots with the common genes
{
    rownames1 <- rownames(df1)
    rownames2 <- rownames(df2)
    common <- intersect(rownames1,rownames2)

    df1 <- df1[common,]
    df2 <- df2[common,]

    return(list(df1,df2))
}
