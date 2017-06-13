GetPs <- function(group_order,n_gene,gene_name,DATA_long)
# This function returns a data frame with the rownames being the group-pair comparison and the column being the gene
{
    #name the groupwise test: i.e. list the two groups being tested
    group_test <- paste(group_order[1],group_order[2],sep=' : ')

    #initialize the data frame
    sig_test_list <- list()
    sig_test_list[[1]] <- data.frame(matrix(nrow=n_gene,ncol=2))
    names(sig_test_list) <- group_test
    rownames(sig_test_list[[1]]) <- gene_name
    colnames(sig_test_list[[1]]) <- c('ratio','p_val')

    SigTestApplyReturn <- apply(gene_name[1:10],1,SigTest)

    #perform the t-tests
    # for (i in 1:n_gene)
    # {
    #     if (i %% 100 == 0){print(i)}
    #     current_gene <- gene_name[i]
    #     indices_group_1 <- DATA_long[,'gene']==gene_name[i] & DATA_long[,'group']==group_order[1]
    #     indices_group_2 <- DATA_long[,'gene']==gene_name[i] & DATA_long[,'group']==group_order[2]
    #     value_group1 <- DATA_long[indices_group_1,'value']
    #     value_group2 <- DATA_long[indices_group_2,'value']
    #     t_test_result <- t.test(value_group1,value_group2,var.equal=TRUE)
    #     p_val <- t_test_result[[3]]
    #     ratio <- median(value_group1)/median(value_group2)
    #     sig_test_list[[1]][current_gene,'ratio'] <- ratio
    #     sig_test_list[[1]][current_gene,'p_val'] <- p_val
    # }

    sig_test_list[[1]][,'log2_ratio'] <- log2(sig_test_list[[1]][,'ratio'])
    sig_test_list[[1]][,'nlog10_p'] <- -log10(sig_test_list[[1]][,'p_val'])

    return(sig_test_list)
}

SigTest(gene_name,groups,DATA_long)
{
    indices_group_1 <- DATA_long[,'gene']==gene_name & DATA_long[,'group']==group_order[1]
    indices_group_2 <- DATA_long[,'gene']==gene_name & DATA_long[,'group']==group_order[2]
    value_group1 <- DATA_long[indices_group_1,'value']
    value_group2 <- DATA_long[indices_group_2,'value']
    t_test_result <- t.test(value_group1,value_group2,var.equal=TRUE)
    p_val <- t_test_result[[3]]
    ratio <- median(value_group1)/median(value_group2)
    SigTestReturn <- list(p_val,ratio)
    return()
}
