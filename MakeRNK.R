MakeRNK <- function(HMR,item)
# Inputs:
#   HMR: is a list of items that is returned from the MakeHeatMap function
#   item: the name of the gene whose correlation with the other genes is used to rank the other genes
#   places the .rnk file one directory above the current working directory in a folder named output

{
    #Extract the corrlelation matrix
    AHI <- HMR[[3]]
    COR <- AHI[[14]]

    #Extract the values of interest and sort them in decreasing order
    values <- sort(COR[,item],decreasing=TRUE)

    #Remove the item correlating with itself
    values <- values[2:length(values)]

    #Place the values in a directory one level above the current directory in a folder called "output"
    output_directory <- StoreHeatmap()
    RNK_file <- paste(output_directory,'RNK_file.rnk',sep='')

    #write the values to the file
    write.table(values,file=RNK_file,sep='\t',row.names=TRUE,quote=FALSE,col.names=FALSE)

return()
}
