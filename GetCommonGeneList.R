GetCommonGeneList <- function(DataDirectories,GeneExpressionFile)
{
    #n_cancers <- length(DataDirectories)
    n_cancers <- 2
    gene_array_list <- list()
    print('getting all gene lists ...')
    for (i in 1:n_cancers)
    {
        print(i)
        GeneExpressionFilePath <- paste(DataDirectories[i],GeneExpressionFile,sep='')
        GeneExpressionInfo <- read.table(file=GeneExpressionFilePath,head=TRUE,sep='\t')
        gene_array_list[[i]] <- GeneExpressionInfo$sample
    }

    print('finding genes in common ...')
    gene_array <- Reduce(intersect,gene_array_list)
    return(gene_array)
}
