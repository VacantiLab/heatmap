assemble_heatmap <- function(AHI)
#This function is needed because heatmap.plus cannot handle only a single grouping scheme, it must be provided at least two
#Thus this function determines if there is one or more than one grouping scheme and uses the appropriate heatmap creating function
#uses dendextend package for coloring dendrogram branches
#uses dynamicTreeCut package to cut the tree for enrichment analysis
#uses enrichR package to calculate cluster enrichments

{
    #unpack the inputs
    GroupColorMatrix <- AHI[[1]]
    DifExpMatx <- AHI[[2]]
    colv <- AHI[[3]]
    rowv <- AHI[[4]]
    break_seq <- AHI[[5]]
    label_rows <- AHI[[6]]
    label_cols <- AHI[[7]]
    output_directory <- AHI[[8]]
    DistanceMethod <- AHI[[9]]
    ClusterMethod <- AHI[[10]]
    C_col <- AHI[[11]]
    C_row <- AHI[[12]]
    Cor_col <- AHI[[13]]
    Cor_row <- AHI[[14]]
    presentation <- AHI[[15]]
    n_clusters <- AHI[[16]]
    graphics_type <- AHI[[17]]
    HeatmapWidth <- AHI[[18]]
    HeatmapHeight <- AHI[[19]]

    library(heatmap.plus)
    n_colors = length(break_seq)-1
    break_seq_0 <- break_seq #save the originally specified break_seq to use in the color key creation
    if (min(DifExpMatx)<break_seq[1]) {break_seq[1]=min(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
    if (max(DifExpMatx)>break_seq[length(break_seq)]) {break_seq[length(break_seq)]=max(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
    heat_map_colors <- colorRampPalette(c('blue','white','red'))(n_colors)
    HeatmapName <- paste(DistanceMethod,'_',ClusterMethod,graphics_type,sep='')
    graphics_file <- paste(output_directory,HeatmapName,sep='')
    graphics_w = HeatmapWidth
    graphics_h = HeatmapHeight

    #color the dendrogram
    color_dend <- FALSE
    if (n_clusters > 1){color_dend = TRUE} # specifies whether or not to color the dendrogram
    cutree_genes = FALSE # initializes the variable containing the cluster membership of genes
    library("RColorBrewer") # loads a library required to select a color palette
    if (color_dend == TRUE)
    {
        if (n_clusters>2){cluster_palette <- brewer.pal(n = n_clusters, name = 'Paired')}
        if (n_clusters==2){cluster_palette <- c('#6B1106','#D18D06')}
        suppressPackageStartupMessages(library(dendextend))
        rowv <- rowv %>% set("branches_k_color", k = n_clusters, value=cluster_palette)
        cutree_genes <- cutree(C_row,n_clusters,order_clusters_as_data=FALSE)
    }

    #open the heatmap graphics file
    if (graphics_type == '.pdf'){pdf(graphics_file,height=graphics_h,width=graphics_w)} #not sure of the units of width and height}
    if (graphics_type == '.png'){png(graphics_file,height=graphics_h,width=graphics_w,units='in',pointsize=24,res=600)}
    if (graphics_type == '.jpeg'){jpeg(graphics_file,height=graphics_h,width=graphics_w,units='in',pointsize=24,res=600)}

    # Redefine parameters passed to heatmap function (the matrix and the column dendrogram) if specified to make a correlation matrix
    if (presentation=='correlation_matrix')
    {
      DifExpMatx <- Cor_row
      colv <- rowv
    }

    if (class(label_rows)=='character')
    {
        genes_to_label <- label_rows
        label_rows <- rownames(DifExpMatx)
        gene_indices <- label_rows %in% genes_to_label
        label_rows[!gene_indices] <- NA
    }

    #make the heat map
    if (!is.null(GroupColorMatrix)){if (dim(GroupColorMatrix)[2]>1)
    {
        heatmap <- heatmap.plus(x=DifExpMatx,
                                Colv=colv,
                                Rowv=rowv,
                                cexRow=0.5,
                                cexCol=0.5,
                                ColSideColors = GroupColorMatrix,
                                col=heat_map_colors,
                                breaks=break_seq,
                                margins=c(6,10), #Can fix the column names being cut off at the bottom
                                scale='none', #This has to be specified to stop the visualization from row scaling (does not impact clustering)
                                labRow=label_rows,
                                labCol=label_cols)
    }}

    if (!is.null(GroupColorMatrix)){if (dim(GroupColorMatrix)[2]==1)
    {
        heatmap <- heatmap(x=DifExpMatx,
                          Colv=colv,
                          Rowv=rowv,
                          cexRow=0.5,
                          cexCol=0.5,
                          ColSideColors = GroupColorMatrix,
                          col=heat_map_colors,
                          breaks=break_seq,
                          #margins=c(5,10), #Can fix the column names being cut off at the bottom
                          scale='none', #This has to be specified to stop the visualization from row scaling (does not impact clustering)
                          labRow=label_rows,
                          labCol=label_cols)
    }}


    if (is.null(GroupColorMatrix))
    {
        library('gplots')
        heatmap <- heatmap.2(x=DifExpMatx,
                          Colv=colv,
                          Rowv=rowv,
                          cexRow=0.5,
                          cexCol=0.5,
                          col=heat_map_colors,
                          breaks=break_seq,
                          #margins=c(5,10), #Can fix the column names being cut off at the bottom
                          scale='none', #This has to be specified to stop the visualization from row scaling (does not impact clustering)
                          labRow=label_rows,
                          labCol=label_cols,
                          trace = 'none',
                          #lhei = c(1,20), #controls the height proportions of the dentrogram with respect to the heat map
                          key = FALSE)
    }

    dev.off() #turn off printing to the specified pdf

    #print the column dendrogram as pdf
    pdf(paste(output_directory,'col_dendrogram',sep=''),height=4,width=10) #not sure of the units of width and height
    plot(C_col,hang=-1)
    dev.off() #turn off printing to the specified pdf


    #print the column dendrogram as pdf
    pdf(paste(output_directory,'row_dendrogram.pdf',sep=''),height=2500,width=300) #not sure of the units of width and height
    plot(rowv,main = "Default colors",lwd=0.5,horiz=TRUE)
    dev.off() #turn off printing to the specified pdf

    #return used variables
    assemble_heatmap_return <- list(heat_map_colors,cutree_genes)
}



# old stuff##########################

# if (color_dend == TRUE)
# {
#     library(dendextend)
#     n_clusters <- 12
#     #par(mfrow = c(1,1))
#     C_row = as.dendrogram(C_row)
#     #requires the dendextend package
#     C_row <- C_row %>% set("branches_k_color", k = n_clusters)
#     #cluster_assignments <- cutree(C_row, k=n_clusters)
#     #cluster_members <- vector(mode='list',length=n_clusters)
#     #enrichment_list <- vector(mode='list',length=n_clusters)
#     #for (i in 1:n_clusters)
#     #{
#     #    current_indices = cluster_assignments == i
#     #    current_members = names(cluster_assignments[current_indices])
#     #    cluster_members[[i]] <- current_members
#     #    current_cluster_members <- cluster_members[[i]]
#     #    current_enrichment <- enrichr(current_cluster_members,'KEGG_2016')
#     #    current_enrichment_df <- as.data.frame(current_enrichment)
#     #    enrichment_list[[i]] <- current_enrichment_df[1:10,]
#
#     #}
# }
#
# dynamic_cut_of_tree = FALSE
# if (dynamic_cut_of_tree == TRUE)
# {
#     dtc <- cutreeDynamic(C_row,minClusterSize=100,deepSplit=0) #requires dynamicTreeCut package
# }
