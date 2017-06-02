GetPs <- function(group_order,n_gene,gene_name,DATA_long)
# This function returns a data frame with the rownames being the group-pair comparison and the column being the gene
{
    #name the groupwise test: i.e. list the two groups being tested
    group_test <- paste(group_order[1],group_order[2],sep=' - ')

    #initialize the data frame
    p_val <- matrix(nrow=1,ncol=n_gene)
    rownames(p_val) <- group_test
    colnames(p_val) <- gene_name

    #perform the t-tests
    for (i in 1:n_gene)
    {
        current_gene <- gene_name[i]
        indices_group_1 <- DATA_long[,'gene']==gene_name[i] & DATA_long[,'group']==group_order[1]
        indices_group_2 <- DATA_long[,'gene']==gene_name[i] & DATA_long[,'group']==group_order[2]
        value_group1 <- DATA_long[indices_group_1,'value']
        value_group2 <- DATA_long[indices_group_2,'value']
        t_test_result <- t.test(value_group1,value_group2,var.equal=TRUE)
        p_val[group_test,current_gene] <- t_test_result[[3]]
    }

    return(p_val)
}
