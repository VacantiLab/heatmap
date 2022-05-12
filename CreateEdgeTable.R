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

    # These are used for filtering purposes
    #     Genes whose correlation with these rows is above a threshold are kept
    coloring_rows = c('PCFlux','UglnM5citrate591','UGlcM3Serine390','UGlcM2Citrate591','UGlnM4Malate419')
    CorrelationThreshold <- 0.85

    #Filter the genes you want to consider: this must be hardcoded in the function below
    FGR <- FilterGenes(DATA,Cor_row,coloring_rows)
    DATA <- FGR[[1]]
    Cor_row <- FGR[[2]]
    gene_names <- FGR[[3]]
    n_genes <- FGR[[4]]

    # Go through the lower left traingle of the correlation matrix and set values equal to TRUE if the meet the correlation threshold
    # Note coloring_row data is manually added to the quantities.txt file if it is not already there (i.e. m3malate or UglnM5citrate591)
    edgeDF <- data.frame(matrix(ncol = 2, nrow = 1))
    colnames(edgeDF) <- c('source','target')
    
    # parameter to set the edge creation criteria
    cor_thresh = CorrelationThreshold
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
    colnames(nodeDF) <- c('ID',coloring_rows[1],'color')

    # Make the node file
    for (i in 1:length(rownames(DATA)))
    {
        print(paste('node: ',i,sep=''))
        current_gene <- rownames(DATA)[i]
        
        # initially have all the nodes be gray
        nodeDF[current_gene,'color'] <- '#FBFBFB'
        nodeDF[current_gene,'ID'] <- current_gene
        nodeDF[current_gene,'z'] <- 1
        
        # specify a color for a correlation with specified rows above threshold 
        colors <- c('#173f5f','#20639b','#3caea3','#B23BFF','#ed553b')
        
        # iterate through each flux to determine if the correlation with the current gene is above threshold
        #     if it is, color the node accordingly
        n_color <- 1
        for (flux in coloring_rows)
        {
          cor_to_gene <- Cor_row[flux,current_gene]
          # parameter to color a gene as correlated with a flux
          if (cor_to_gene >= CorrelationThreshold)
          {nodeDF[current_gene,'color'] <- colors[n_color]
           nodeDF[current_gene,'z'] <- 2
          }
          n_color <- n_color+1
        }
        
    }
    # put the higher z-score nodes at the end of the data frame
    #     Gephi will draw the nodes in the order in which they appear
    #     If you want the higher z-scores on top, they need to be drawn last
    #     The value of z is not interpreted byu Gephi (at least right now)
    nodeDF <- nodeDF[order(nodeDF$z),]
    
    # Remove the nodes and edges containing the metabolic fluxes
    #     This is a protein-only network, the fluxes are present only for correlation calculations
    NodeRowsToKeep <- !(rownames(nodeDF) %in% coloring_rows)
    nodeDF <- nodeDF[NodeRowsToKeep,]
    
    EdgeRowsToKeep <- !((edgeDF[,'source'] %in% coloring_rows) | (edgeDF[,'target'] %in% coloring_rows))
    edgeDF <- edgeDF[EdgeRowsToKeep,]
    
    #write the tables required to make the graph file (gdf file) made by make_gdf.py
    write_directory <- '/Users/nate/Desktop/temporary/'
    write.table(edgeDF,file=paste(write_directory,'edges.csv',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")
    write.table(nodeDF[,c('ID','color')],file=paste(write_directory,'nodes.csv',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")
}




FilterGenes <- function(DATA,Cor_row,coloring_rows)
{

  DATA_matrix <- as.matrix(DATA)
  gene_names <- rownames(Cor_row)
  n_genes <- length(gene_names)
  samples <- colnames(DATA_matrix)


  #Filter genes
  selection_criteria <- rep(0,n_genes)
  for (i in 1:n_genes)
  {
     # Set the differential expression threshold
      DifferenceThreshold <- 0.1
      selection_criteria[i] <- (sum(DATA[i,] >= DifferenceThreshold) >= 2) | (gene_names[i] %in% coloring_rows)

      # If specified, keep genes whose correlation is above a threshold with the color indicator gene/MID value even if they do not meet the differential expression threshold
      #     The threshold needs to be met in a correlation with one of the specified fluxes in coloring_rows
      AdditionalCorrelatingGenesCriteria <- FALSE
      CorFilterCutOff <- 1
      CriteriaVector <- abs(Cor_row[coloring_rows,i]) >= CorFilterCutOff
      if (sum(CriteriaVector) >= 1 & AdditionalCorrelatingGenesCriteria)
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
