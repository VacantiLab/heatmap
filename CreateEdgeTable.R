CreateEdgeTable <- function(Cor_row,DATA)
{
    select_samples <- c('bt20_rotenone','hcc1419_rotenone','mcf7_rotenone','t47d_rotenone','mdamb157_rotenone')
    DATA <- DATA[,select_samples]
    DATA_matrix <- as.matrix(DATA)
    gene_names <- rownames(Cor_row)
    n_genes <- length(gene_names)
    rel_var_array <- rep(0,n_genes)
    for (i in 1:n_genes)
    {
        rel_var_array[i] <- var(DATA_matrix[i,])/mean(abs(DATA_matrix[i,]))
    }
    rel_var_thresh <- 0
    rel_var_thresh_logical <- rel_var_array > rel_var_thresh

    gene_names <- gene_names[rel_var_thresh_logical]
    Cor_row <- Cor_row[gene_names,gene_names]
    n_genes <- length(gene_names)

    EdgeDF <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(EdgeDF) <- c('source','target')
    cor_thresh = 0.8
    Cor_thresh_logical <- Cor_row
    Cor_thresh_logical[,] <- FALSE
    for (i in 2:n_genes)
    {
        print(paste('gene: ',i,sep=''))
        for (j in (1:(i-1)))
        {
            if (Cor_row[i,j] > cor_thresh)
            {
                Cor_thresh_logical[i,j] <- TRUE
            }
        }
    }

    edge_iterator <- 1
    n_interactions <- sum(Cor_thresh_logical)
    EdgeDF <- data.frame(matrix(0,n_interactions,2))
    colnames(EdgeDF) <- c('source','target')
    edge_iterator <- 1
    for (i in 2:n_genes)
    {
        print(paste('gene: ',i,sep=''))
        current_gene_mates <- gene_names[Cor_thresh_logical[i,]==1]
        if (length(current_gene_mates > 0))
        {
            EdgeDF[edge_iterator:(edge_iterator+length(current_gene_mates)-1),1] = rep(gene_names[i],length(current_gene_mates))
            EdgeDF[edge_iterator:(edge_iterator+length(current_gene_mates)-1),2] = current_gene_mates
            edge_iterator <- edge_iterator + length(current_gene_mates)
        }

        # for (j in (1:(i-1)))
        # {
        #     current_gene_mates <- gene_names[Cor_thresh_logical[]]
        #     if (Cor_thresh_logical[i,j] == 1)
        #     {
        #         EdgeDF[edge_iterator,'source'] <- gene_names[i]
        #         EdgeDF[edge_iterator,'target'] <- gene_names[j]
        #         edge_iterator <- edge_iterator + 1
        #     }
        # }
    }

    browser()
}
