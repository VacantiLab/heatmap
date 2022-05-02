CreateEdgeTable <- function(Cor_row,DATA)
# This function is run by MakeHeatMap when the visualization is set to 'edge_table'
#     It creates the edge and node files needed to make a .gdf file readable by Gephi
#     The .gdf file is made using the node and edge files and is assembled by a separate python funcion in the heatap package
# The following are hard-coded parameters:
#     coloring_row: sets the row of DATA which the color of the points in the network are based off
#     DifferenceThreshold: sets the differential threshold of treatment vs. no treatment for a point to be accepted
#         This assumes treated vs. not treated are next to each other and there are 10 samples
#     CorFilterCutOff: This establishes a correlation magnitude with coloring_row above which rows are not subject to the DifferenceThreshold

{

    # The gene/row whose correlation relationship with the others sets the color of the network
    #     This is hard-coded for now
    coloring_row = 'UglnM5citrate591'

    #Filter the genes you want to consider: this must be hardcoded in the function below
    FGR <- FilterGenes(DATA,Cor_row,coloring_row)
    DATA <- FGR[[1]]
    Cor_row <- FGR[[2]]
    gene_names <- FGR[[3]]
    n_genes <- FGR[[4]]

    # Go through the lower left traingle of the correlation matrix and set values equal to TRUE if the meet the correlation threshold
    # Note coloring_row data is manually added to the quantities.txt file if it is not already there (i.e. m3malate or UglnM5citrate591)
    edgeDF <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(edgeDF) <- c('source','target')
    cor_thresh = 0.875
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
    colnames(nodeDF) <- c('ID',coloring_row,'color')
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
        cor_to_m3mal <- Cor_row[coloring_row,current_gene]
        nodeDF[current_gene,coloring_row] <- cor_to_m3mal
        fraction_up_scale <- (cor_to_m3mal - (-1))/2
        color_index <- round(1 + fraction_up_scale*(n_colors-1))
        nodeDF[current_gene,'color'] <- custom_palette[color_index]
    }

    #write the tables required to make the graph file (gdf file) made by make_gdf.py
    write_directory <- '/Users/nate/Desktop/temporary/'
    write.table(edgeDF,file=paste(write_directory,'edges.csv',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")
    write.table(nodeDF[,c('ID','color')],file=paste(write_directory,'nodes.csv',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")
}

FilterGenes <- function(DATA,Cor_row,coloring_row)
{

  DATA_matrix <- as.matrix(DATA)
  gene_names <- rownames(Cor_row)
  n_genes <- length(gene_names)
  samples <- colnames(DATA_matrix)


  #Filter genes
  selection_criteria <- rep(0,n_genes)
  for (i in 1:n_genes)
  {
      Difference1 <- abs(DATA_matrix[i,samples[1]]-DATA_matrix[i,samples[2]])
      Difference2 <- abs(DATA_matrix[i,samples[3]]-DATA_matrix[i,samples[4]])
      Difference3 <- abs(DATA_matrix[i,samples[5]]-DATA_matrix[i,samples[6]])
      Difference4 <- abs(DATA_matrix[i,samples[7]]-DATA_matrix[i,samples[8]])
      Difference5 <- abs(DATA_matrix[i,samples[9]]-DATA_matrix[i,samples[10]])

      # Set the differential expression threshold
      DifferenceThreshold <- 0.65
      selection_criteria[i] <- (Difference1 >= DifferenceThreshold) | (Difference2 >= DifferenceThreshold) | (Difference3 >= DifferenceThreshold) | (Difference4 >= DifferenceThreshold) | (Difference5 >= DifferenceThreshold)

      # Keep genes whose correlation is above a threshold with the color indicator gene/MID value even if they do not meet the differential expression threshold
      CorFilterCutOff <- 0.93
      if (abs(Cor_row[coloring_row,i]) > CorFilterCutOff)
      {
          selection_criteria[i] <- TRUE
      }
  }
  selection_criteria <- as.logical(selection_criteria)

  gene_names <- gene_names[selection_criteria]
  Cor_row <- Cor_row[gene_names,gene_names]
  DATA <- DATA[gene_names,]
  n_genes <- length(gene_names)

  FGR <- list(DATA,
              Cor_row,
              gene_names,
              n_genes)

  return(FGR)

}
