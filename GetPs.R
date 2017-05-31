GetPs <- function(group_order_original,n_gene,gene_name,DATA_long)
{
    group_test <- paste(group_order_original[1],group_order_original[2],sep=' - ')
    p_val_df <- data.frame(matrix(nrow=1,ncol=n_gene))
    rownames(p_val_df) <- group_test
    colnames(p_val_df) <- gene_name
    for (i in 1:n_gene)
    {
        current_gene <- gene_name[i]
        indices_group_1 <- DATA_long[,'gene']==gene_name[i] & DATA_long[,'group']==group_order_original[1]
        indices_group_2 <- DATA_long[,'gene']==gene_name[i] & DATA_long[,'group']==group_order_original[2]
        value_group1 <- DATA_long[indices_group_1,'value']
        value_group2 <- DATA_long[indices_group_2,'value']
        t_test_result <- t.test(value_group1,value_group2,var.equal=TRUE)
        p_val_df[group_test,current_gene] <- t_test_result[3]
    }

    return(p_val_df)
}
