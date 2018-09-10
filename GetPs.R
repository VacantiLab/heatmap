GetPs <- function(group_order,gene_name,DATA,groups_corresponding,visualization)
# This function returns a data frame with the rownames being the group-pair comparison and the column being the gene
{
    #name the groupwise test: i.e. list the two groups being tested
    group_test <- paste(group_order[1],group_order[2],sep=' : ')

    #remove genes from the list to consdier if the data does not exist (a t-test will not be attempted)
    list_to_consider <- gene_name %in% rownames(DATA)
    gene_name <- gene_name[list_to_consider]

    #initialize the data frame
    sig_test_list <- list()
    n_genes <- length(gene_name)
    sig_test_list[[1]] <- data.frame(matrix(nrow=n_genes,ncol=2))
    names(sig_test_list) <- group_test
    rownames(sig_test_list[[1]]) <- gene_name
    colnames(sig_test_list[[1]]) <- c('ratio','p_val')

    gene_name <- as.matrix(gene_name,nrow=n_genes)

    sig_test_list <- list()
    sig_test_list[[1]] <- data.frame(matrix(nrow=n_genes,ncol=2))
    names(sig_test_list) <- group_test
    rownames(sig_test_list[[1]]) <- gene_name[,1]
    colnames(sig_test_list[[1]]) <- c('ratio','p_val')

    group1_members <- rownames(groups_corresponding[groups_corresponding[,1]==group_order[1],1,drop=FALSE])
    group2_members <- rownames(groups_corresponding[groups_corresponding[,1]==group_order[2],1,drop=FALSE])
    SigTestApplyReturn <- t(apply(gene_name,1,SigTest,DATA=DATA,group1_members=group1_members,group2_members=group2_members,visualization=visualization))
    rownames(SigTestApplyReturn) <- gene_name
    colnames(SigTestApplyReturn) <- c('p_val','ratio')

    rownames(SigTestApplyReturn) <- gene_name
    colnames(SigTestApplyReturn) <- c('p_val','ratio')

    sig_test_list[[1]][gene_name,'p_val'] <- SigTestApplyReturn[gene_name,'p_val']
    sig_test_list[[1]][gene_name,'ratio'] <- SigTestApplyReturn[gene_name,'ratio']
    sig_test_list[[1]][,'log2_ratio'] <- log2(sig_test_list[[1]][,'ratio']) #this will produce NA and a warning if the ratio is negative (could be depending on the transformation)
    sig_test_list[[1]][,'nlog10_p'] <- -log10(sig_test_list[[1]][,'p_val'])
    sig_test_list[[1]] <- sig_test_list[[1]][c('p_val','nlog10_p','ratio','log2_ratio')] #reorders data frame by column name

    return(sig_test_list)
}

#Supporting function to calculate p-values
SigTest <- function(gene_name,DATA,group1_members,group2_members,visualization)
{
    value_group1 <- as.numeric(DATA[gene_name,group1_members])
    value_group2 <- as.numeric(DATA[gene_name,group2_members])
    #alternative_option = 'greater'
    #alternative_option = 'less'
    alternative_option = 'two.sided'

    # perform the t-tests on log2 transformed values if specified
    #     this true for volcano plots and not boxplots!
    #     volcano plots never have a transformation specified becuase that is done automatically in the program
    value_group1_transformed = value_group1
    value_group2_transformed = value_group2
    if (visualization == 'volcanoplot')
    {
        value_group1_transformed = log2(value_group1)
        value_group2_transformed = log2(value_group2)
    }

    t_test_result <- t.test(value_group1_transformed,value_group2_transformed,var.equal=TRUE,alternative=alternative_option)
    p_val <- t_test_result[[3]]
    ratio <- median(value_group1)/median(value_group2)
    SigTestReturn <- c(p_val,ratio)
    return(SigTestReturn)
}
