assemble_heatmap <- function(GroupColorMatrix,DifExpMatx,colv,rowv,heat_map_colors,break_seq,label_rows)
#This function is needed because heatmap.plus cannot handle only a single grouping scheme, it must be provided at least two
#Thus this function determines if there is one or more than one grouping scheme and uses the appropriate heatmap creating function
  
{
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
                            labRow=label_rows)
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
                      labRow=label_rows)
  }
  
}