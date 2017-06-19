assemble_heatmap <- function(GroupColorMatrix,DifExpMatx,colv,rowv,break_seq,label_rows,label_cols,HeatmapDirectory,DistanceMethod,ClusterMethod)
#This function is needed because heatmap.plus cannot handle only a single grouping scheme, it must be provided at least two
#Thus this function determines if there is one or more than one grouping scheme and uses the appropriate heatmap creating function

{
    n_colors = length(break_seq)-1
    break_seq_0 <- break_seq #save the originally specified break_seq to use in the color key creation
    if (min(DifExpMatx)<break_seq[1]) {break_seq[1]=min(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
    if (max(DifExpMatx)>break_seq[length(break_seq)]) {break_seq[length(break_seq)]=max(DifExpMatx)} #This needs to be done because heatmap.plus assigns white to everything outside the range
    heat_map_colors <- colorRampPalette(c('blue','white','red'))(n_colors)
    HeatmapName <- paste(DistanceMethod,'_',ClusterMethod,'.pdf',sep='')
    PDF_file <- paste(HeatmapDirectory,HeatmapName,sep='')
    PdfW = 7
    PdfH = 7
    pdf(PDF_file,height=PdfH,width=PdfW) #not sure of the units of width and height

    if (dim(GroupColorMatrix)[2]>1)
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
    }

    if (dim(GroupColorMatrix)[2]==1)
    {
        heatmap <- heatmap(x=DifExpMatx,
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
    }

    #return used variables
    assemble_heatmap_return <- list(heat_map_colors)
}
