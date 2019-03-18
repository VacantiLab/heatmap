CreateEdgeTable <- function(Cor_row,DATA)
{
    #Filter the genes you want to consider: this must be hardcoded
    select_samples <- c('bt20_rotenone','hcc1419_rotenone','mcf7_rotenone','t47d_rotenone','mdamb157_rotenone')
    DATA <- DATA[,select_samples] #this is only used to select genes. All columns are used in the correlation filter (already included in Cor_row)
    DATA_matrix <- as.matrix(DATA)
    gene_names <- rownames(Cor_row)
    n_genes <- length(gene_names)

    #Filter genes
    selection_criteria <- rep(0,n_genes)
    for (i in 1:n_genes)
    {
        selection_criteria[i] <- (max(DATA_matrix[i,]) > 0.2) & (min(DATA_matrix[i,]) < -0.2)
    }
    selection_criteria <- as.logical(selection_criteria)

    gene_names <- gene_names[selection_criteria]
    Cor_row <- Cor_row[gene_names,gene_names]
    DATA <- DATA[gene_names,]
    n_genes <- length(gene_names)

    # Go through the lower left traingle of the correlation matrix and set values equal to TRUE if the meet the correlation threshold
    # Note m3malate data is hardcoded into the DATA matrix by the ArrangeData function if visualization is edge_table
    edgeDF <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(edgeDF) <- c('source','target')
    cor_thresh = 0.85
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
    edgeDF <- data.frame(matrix(0,n_interactions,2))
    colnames(edgeDF) <- c('source','target')
    edge_iterator <- 1
    for (i in 2:n_genes)
    {
        print(paste('record_edge: ',i,sep=''))
        current_gene_mates <- gene_names[Cor_thresh_logical[i,]==1]
        if (length(current_gene_mates > 0))
        {
            edgeDF[edge_iterator:(edge_iterator+length(current_gene_mates)-1),1] = rep(gene_names[i],length(current_gene_mates))
            edgeDF[edge_iterator:(edge_iterator+length(current_gene_mates)-1),2] = current_gene_mates
            edge_iterator <- edge_iterator + length(current_gene_mates)
        }
    }

    nodeDF <- data.frame(matrix(0,nrow(DATA),3))
    rownames(nodeDF) <- rownames(DATA)
    colnames(nodeDF) <- c('ID','m3malcor','color')
    colfunc <- colorRampPalette(c('#000055',"#006CFF","white","#FF6B6B",'#4A0000'))
    n_colors <- 20
    custom_palette <- colfunc(n_colors)
    custom_palette <- c(custom_palette)
    n_colors <- n_colors
    # Make the node file
    for (i in 1:length(rownames(DATA)))
    {
        print(paste('node: ',i,sep=''))
        current_gene <- rownames(DATA)[i]
        nodeDF[current_gene,'ID'] <- current_gene
        cor_to_m3mal <- Cor_row['m3malate',current_gene]
        nodeDF[current_gene,'m3malcor'] <- cor_to_m3mal
        fraction_up_scale <- (cor_to_m3mal - (-1))/2
        color_index <- round(1 + fraction_up_scale*(n_colors-1))
        nodeDF[current_gene,'color'] <- custom_palette[color_index]
    }

    #write the tables required to make the graph file (gdf file) made by make_gdf.py
    write_directory <- '/Users/nate/Desktop/temporary/'
    write.table(edgeDF,file=paste(write_directory,'edges.csv',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")
    write.table(nodeDF[,c('ID','color')],file=paste(write_directory,'nodes.csv',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")

    browser()
}
