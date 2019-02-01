assemble_heatmap <- function(GroupColorMatrix,DifExpMatx,colv,rowv,break_seq,label_rows,label_cols,HeatmapDirectory,DistanceMethod,ClusterMethod,C_col,C_row,Cor_col,Cor_row,presentation)
#This function is needed because heatmap.plus cannot handle only a single grouping scheme, it must be provided at least two
#Thus this function determines if there is one or more than one grouping scheme and uses the appropriate heatmap creating function
#uses dendextend package for coloring dendrogram branches
#uses dynamicTreeCut package to cut the tree for enrichment analysis
#uses enrichR package to calculate cluster enrichments

{
    library(heatmap.plus)
    n_colors = length(break_seq)-1
    break_seq_0 <- break_seq #save the originally specified break_seq to use in the color key creation
    if (min(DifExpMatx)<break_seq[1]) {break_seq[1]=min(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
    if (max(DifExpMatx)>break_seq[length(break_seq)]) {break_seq[length(break_seq)]=max(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
    heat_map_colors <- colorRampPalette(c('blue','white','red'))(n_colors)
    graphics_type <- '.pdf'
    HeatmapName <- paste(DistanceMethod,'_',ClusterMethod,graphics_type,sep='')
    graphics_file <- paste(HeatmapDirectory,HeatmapName,sep='')
    graphics_w = 8
    graphics_h = 8

    #color the dendrogram
    color_dend = TRUE
    if (color_dend == TRUE)
    {
        library(dendextend)
        n_clusters <- 12
        rowv <- rowv %>% set("branches_k_color", k = n_clusters)
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
                                margins=c(5,10),
                                scale='none', #This has to be specified in heatmap.plus and heatmap but not in heatmap.2
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
                          #margins=c(5,10),
                          scale='none', #This has to be specified in heatmap.plus and heatmap but not in heatmap.2
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
                          #margins=c(5,10),
                          scale='none', #This has to be specified in heatmap.plus and heatmap but not in heatmap.2
                          labRow=label_rows,
                          labCol=label_cols,
                          trace = 'none',
                          #lhei = c(1,20), #controls the height proportions of the dentrogram with respect to the heat map
                          key = FALSE)
    }

    dev.off() #turn off printing to the specified pdf

    #print the column dendrogram as pdf
    pdf(paste(HeatmapDirectory,'col_dendrogram',sep=''),height=4,width=10) #not sure of the units of width and height
    plot(C_col,hang=-1)
    dev.off() #turn off printing to the specified pdf

    #print the column dendrogram as pdf
    pdf(paste(HeatmapDirectory,'row_dendrogram.pdf',sep=''),height=2500,width=100) #not sure of the units of width and height
    plot(rowv,main = "Default colors",lwd=0.5,horiz=TRUE)
    dev.off() #turn off printing to the specified pdf

    #return used variables
    assemble_heatmap_return <- list(heat_map_colors)
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
