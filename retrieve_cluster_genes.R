retrieve_cluster_genes <- function(ctg,cluster_id)
{
    #retrieves the genes which are members of the specified cluster
    #cg is a vector output from cutree(), an output from MakeHeatmap
    #cluster ID is an integer representing the desired cluster
    #    the cluster ID can be obtained from the cluster_color_legend.pdf
    cluster_genes <- names(ctg[ctg==cluster_id])
    return(cluster_genes)
}
