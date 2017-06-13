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

    n_genes <- length(gene_name)
    gene_name <- as.matrix(gene_name,nrow=n_genes)

    sig_test_list <- list()
    sig_test_list[[1]] <- data.frame(matrix(nrow=n_genes,ncol=2))
    names(sig_test_list) <- group_test
    rownames(sig_test_list[[1]]) <- gene_name[,1]
    colnames(sig_test_list[[1]]) <- c('ratio','p_val')

    start_time <- Sys.time()
    SigTestApplyReturn <- t(apply(gene_name,1,SigTest,group_order=group_order,DATA_long=DATA_long))
    end_time <- Sys.time()
    print(paste('p_value and ratio calculation time:',' ',end_time-start_time,' ','minutes',sep=''))

    rownames(SigTestApplyReturn) <- gene_name
    colnames(SigTestApplyReturn) <- c('p_val','ratio')

    sig_test_list[[1]][gene_name,'p_val'] <- SigTestApplyReturn[gene_name,'p_val']
    sig_test_list[[1]][gene_name,'ratio'] <- SigTestApplyReturn[gene_name,'ratio']
    sig_test_list[[1]][,'log2_ratio'] <- log2(sig_test_list[[1]][,'ratio'])
    sig_test_list[[1]][,'nlog10_p'] <- -log10(sig_test_list[[1]][,'p_val'])

    return(sig_test_list)
}

SigTest <- function(gene_name,group_order,DATA_long)
{
    value_group1 <- with(DATA_long,value[group==group_order[1] & gene==gene_name])
    value_group2 <- with(DATA_long,value[group==group_order[2] & gene==gene_name])
    t_test_result <- t.test(value_group1,value_group2,var.equal=TRUE)
    p_val <- t_test_result[[3]]
    ratio <- median(value_group1)/median(value_group2)
    SigTestReturn <- c(p_val,ratio)
    return(SigTestReturn)
}
