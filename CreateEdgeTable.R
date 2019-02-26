CreateEdgeTable <- function(Cor_row)
{
    EdgeDF <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(EdgeDF) <- c('source','target')
    cor_thresh = 0.5
    Cor_thresh_logical <- Cor_row
    Cor_thresh_logical[,] <- FALSE
    gene_names <- rownames(Cor_row)
    n_genes <- length(gene_names)
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
    for (i in 2:n_genes)
    {
        print(paste('gene: ',i,sep=''))
        current_gene_mates <- gene_names[Cor_thresh_logical[i,]==1]
        if (length(current_gene_mates > 0))
        {
            EdgeDF[1:length(current_gene_mates),1] = rep(gene_names[i],length(current_gene_mates))
            EdgeDF[1:length(current_gene_mates),2] = current_gene_mates
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
