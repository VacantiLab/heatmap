probes_to_genes <- function(probe_txt)

{
    #load the propper library
    library('hgu133plus2.db')

    #read the data from the file into a data frame using a user-defined function
    DATA <- read_txt_to_df(probe_txt)

    #the rownames of the data frame are the microarray probes
    ids <- rownames(DATA)

    #get a key for the probes using the specific microarray database from bioclite
    #    note bioclite must be installed along with the library for the plate
    #    plates are standardized and there libraries can be installed
    key1 <- select(hgu133plus2.db,ids,c('SYMBOL','GENENAME'),'PROBEID')
    counter <- 0
    for (i in rownames(DATA))
    {
        #if the probe has a gene matched to it, record that gene
        if (i %in% key1[,'PROBEID'])
        {
            loc <- match(i,key1[,'PROBEID'])
            DATA[i,'gene_symbol'] <- key1[loc,'SYMBOL']
        }

        #if the probe does not match to a gene, get rid of the corresponding row
        #    some probes are internal standards for the assay
        if (!(i %in% key1[,'PROBEID']))
        {
            DATA <- DATA[!rownames(DATA) %in% i,]
        }

        #print progress on seach for gene symbols
        counter <- counter + 1
        if (counter %% 1000 == 0)
        {
            print(paste('fining gene corresponding to probe number: ', toString(counter), sep=''))
        }
    }

    #average the row values of probes that mapped to identical gene symbols
    #    the gene symbols are stored as 'Group.1'
    #    the 'gene_symbol' column will be converted to NA because it contains strings
    #        this produces warnings that can be ignored
    DATA <- aggregate(DATA,by=list(DATA[,'gene_symbol']),FUN=mean,na.rm=TRUE)
    DATA[,'gene_symbol'] <- NULL
    rownames(DATA) <- DATA[,'Group.1']
    DATA[,'Group.1'] <- NULL

    #the returned data frame has gene symbols as row names and sample names as column names
    return(DATA)
}

###############################################################################

read_txt_to_df <- function(txt_directory)
#Function to read a tab delimited text file into a data frame
{
    DF <- read.table(file=txt_directory,head=TRUE,check.names=FALSE,sep='\t') #check.names=FALSE prevents an 'X' from being added to the numeric column names
    #Name the rows of the data frame as the genes given in the first column of the data frame
    RowNames <- as.character(DF[,1])
    rownames(DF) <- RowNames
    DF[,1] <- NULL #remove the first column of the data frame as it is no longer needed
    return(DF)
}

###############################################################################
