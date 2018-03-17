assemble_heatmap <- function(GroupColorMatrix,DifExpMatx,colv,rowv,break_seq,label_rows,label_cols,HeatmapDirectory,DistanceMethod,ClusterMethod,C_col,C_row)
#This function is needed because heatmap.plus cannot handle only a single grouping scheme, it must be provided at least two
#Thus this function determines if there is one or more than one grouping scheme and uses the appropriate heatmap creating function
#uses dendextend package for coloring dendrogram branches
#uses dynamicTreeCut package to cut the tree for enrichment analysis

{
    n_colors = length(break_seq)-1
    break_seq_0 <- break_seq #save the originally specified break_seq to use in the color key creation
    if (min(DifExpMatx)<break_seq[1]) {break_seq[1]=min(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
    if (max(DifExpMatx)>break_seq[length(break_seq)]) {break_seq[length(break_seq)]=max(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
    heat_map_colors <- colorRampPalette(c('blue','white','red'))(n_colors)
    graphics_type <- '.png'
    HeatmapName <- paste(DistanceMethod,'_',ClusterMethod,graphics_type,sep='')
    graphics_file <- paste(HeatmapDirectory,HeatmapName,sep='')
    graphics_w = 1000
    graphics_h = 1000

    if (graphics_type == '.pdf')
    {
        pdf(graphics_file,height=4,width=4) #not sure of the units of width and height
    }

    if (graphics_type == '.png')
    {
        png(graphics_file,height=graphics_h,width=graphics_w,units='px',pointsize=24)
    }

    if (!is.null(GroupColorMatrix)){if (dim(GroupColorMatrix)[2]>1)
    {
        heatmap <- heatmap.plus(x=DifExpMatx,
                                Colv=colv,
                                Rowv=rowv,
                                cexRow=1,
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
                          margins=c(5,10),
                          scale='none', #This has to be specified in heatmap.plus and heatmap but not in heatmap.2
                          labRow=label_rows,
                          labCol=label_cols)
    }}

    if (is.null(GroupColorMatrix))
    {
        heatmap <- heatmap(x=DifExpMatx,
                          Colv=colv,
                          Rowv=rowv,
                          cexRow=0.5,
                          cexCol=0.5,
                          col=heat_map_colors,
                          breaks=break_seq,
                          margins=c(5,10),
                          scale='none', #This has to be specified in heatmap.plus and heatmap but not in heatmap.2
                          labRow=label_rows,
                          labCol=label_cols)
    }

    dev.off() #turn off printing to the specified pdf

    #print the column dendrogram as pdf
    pdf(paste(HeatmapDirectory,'col_dendrogram',sep=''),height=4,width=4) #not sure of the units of width and height
    plot(C_col,hang=-1)
    dev.off() #turn off printing to the specified pdf

    #print the column dendrogram as pdf
    color_dend = FALSE
    pdf(paste(HeatmapDirectory,'row_dendrogram',sep=''),height=10,width=10) #not sure of the units of width and height
    if (color_dend == TRUE)
    {
        par(mfrow = c(1,1))
        C_row = as.dendrogram(C_row)
        #requires the dendextend package
        C_row %>% set("branches_k_color", k = 10) %>% plot(main = "Row Dendrogram")
    }

    dynamic_cut_of_tree = FALSE
    if (dynamic_cut_of_tree == TRUE)
    {
        dtc <- cutreeDynamic(C_row,minClusterSize=100,deepSplit=0) #requires dynamicTreeCut package
    }

    if (color_dend == FALSE)
    {
        plot(C_row,hang=-1,lwd=0.5)
    }
    dev.off() #turn off printing to the specified pdf

    #return used variables
    assemble_heatmap_return <- list(heat_map_colors)
}
