CreateEdgeTable <- function(Cor_row,DATA)
{
    #Filter the genes you want to consider: this must be hardcoded
    select_samples <- c('bt20_rotenone','hcc1419_rotenone','mcf7_rotenone','t47d_rotenone','mdamb157_rotenone')
    DATA <- DATA[,select_samples]
    DATA_matrix <- as.matrix(DATA)
    gene_names <- rownames(Cor_row)
    n_genes <- length(gene_names)

    DATA['m3malate',c('bt20_rotenone','hcc1419_rotenone','mcf7_rotenone','t47d_rotenone','mdamb157_rotenone')] <- log2(c(1.87,0.575,0.515,0.300,1.26))

    # selection_criteria <- rep(0,n_genes)
    # for (i in 1:n_genes)
    # {
    #     selection_criteria[i] <- max(abs(DATA_matrix[i,])) ==  abs(DATA_matrix[i,'bt20_rotenone'])
    # }
    # selection_criteria <- as.logical(selection_criteria)
    #
    # gene_names <- gene_names[selection_criteria]
    # Cor_row <- Cor_row[gene_names,gene_names]
    # n_genes <- length(gene_names)


    # Go through the lower left traingle of the correlation matrix and set values equal to TRUE if the meet the correlation threshold
    EdgeDF <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(EdgeDF) <- c('source','target')
    cor_thresh = 0.95
    Cor_thresh_logical <- Cor_row
    Cor_thresh_logical[,] <- FALSE
    for (i in 2:n_genes)
    {
        print(paste('check_edge: ',i,sep=''))
        for (j in (1:(i-1)))
        {
            if (Cor_row[i,j] > cor_thresh)
            {
                Cor_thresh_logical[i,j] <- TRUE
            }
        }
    }

    # Make the edge data frame
    #     For all of the TRUEs in the above matrix, record the associated gene pair in a row which represents the edge
    edge_iterator <- 1
    n_interactions <- sum(Cor_thresh_logical)
    EdgeDF <- data.frame(matrix(0,n_interactions,2))
    colnames(EdgeDF) <- c('source','target')
    edge_iterator <- 1
    for (i in 2:n_genes)
    {
        print(paste('record_edge: ',i,sep=''))
        current_gene_mates <- gene_names[Cor_thresh_logical[i,]==1]
        if (length(current_gene_mates > 0))
        {
            EdgeDF[edge_iterator:(edge_iterator+length(current_gene_mates)-1),1] = rep(gene_names[i],length(current_gene_mates))
            EdgeDF[edge_iterator:(edge_iterator+length(current_gene_mates)-1),2] = current_gene_mates
            edge_iterator <- edge_iterator + length(current_gene_mates)
        }
    }

    nodesDF <- data.frame(matrix(0,nrow(DATA),3))
    rownames(nodesDF) <- rownames(DATA)
    colnames(nodesDF) <- c('ID','m3malcor','color')
    colfunc <- colorRampPalette(c("blue", "red"))
    custom_palette <- colfunc(11)
    # Make the node file
    for (i in 1:length(rownames(DATA)))
    {
        print(paste('node: ',i,sep=''))
        current_gene <- rownames(DATA)[i]
        nodesDF[current_gene,'ID'] <- current_gene
        cor_to_m3mal <- cor(as.matrix(DATA)[current_gene,],as.matrix(DATA)['m3malate',],method='pearson')
        nodesDF[current_gene,'m3malcor'] <- cor_to_m3mal
        color_index <- round(cor_to_m3mal*5+6)
        nodesDF[current_gene,'color'] <- custom_palette[color_index]
    }

    browser()
}
