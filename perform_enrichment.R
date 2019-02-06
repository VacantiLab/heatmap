perform_enrichment <- function(hm,cluster_id,list)
{
    ctg <- hm[[5]]
    if (!('enrichR' %in% .packages())){library(enrichR)}
    cluster_genes <- names(ctg[ctg==cluster_id])
    enrichment_output <- enrichr(cluster_genes,list)
    enrichment_output <- as.data.frame(enrichment_output)
    df1 <- enrichment_output[,c(paste(list,'.Term',sep=''),paste(list,'.Overlap',sep=''),paste(list,'.P.value',sep=''),paste(list,'.Adjusted.P.value',sep=''))][1:20,]
    df2 <- enrichment_output[,paste(list,'.Genes',sep=''),drop=FALSE][1:20,,drop=FALSE]
    enrichment_output <- list(df1,df2)
    return(enrichment_output)
}
