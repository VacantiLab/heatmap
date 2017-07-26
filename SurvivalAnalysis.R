SurvivalAnalysis <- function()
# This function is currently being modified
# This function analyzes patient survival and gene expression data over different tumor types
# The data contained for each tumor type are stored in separate directories from the other tumor types
#     Each tumor-type directory contains a text file with gene expression data and separate text file with patient survival data
# Function design:
#     This function loops through the tumor-types by looping through the tumor-type directories
#     Once in a tumor-type directory, the function loops through the genes whose expressions were measured in the tumor-type study
#     For each gene, a two parameter Cox survival model is fit, where patient age and expression of the gene are the Cox parameters
#     For each gene within the tumor-type, a hazard coefficient along with it's 95% CI are calculated
#     If the hazard coefficient is above 1 without the 95% CI overlapping 1, the gene is considered a risk factor and the variable risk_factor is given a value of 1
#     If the hazard coefficient is below 1 without the 95% CI overlapping 1, the gene is considered a benefitial and the variable risk_factor is given a value of -1
#     If the 95% CI overlaps 1, the gene is considered tumor neutral and the variable risk_factor is given a value of 0
#     Once the nested loops through tumor-type and genes are completed, each gene for each tumor-type will have a risk_factor variable associated with it
#     This will be stored in a data frame
#          Euclidian distance hierarchical clustering of the genes will allow visualization of which genes act in concert as risk or survival factors
#          Similar clustering of the tumors will allow visualization of which tumors share common risk and survival factor profiles
{
    #load requisite libraries
    library(ggplot2)
    library(survival)
    library(tidyr)
    library(matrixStats)

    CancerIdentifiers <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','COADREAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC',
                           'LUAD','LUNG','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')

    #KICH moved to first because it was having problems
    #CancerIdentifiers <- c('KICH','ACC','BLCA','BRCA','CESC','CHOL','COAD','COADREAD','DLBC','ESCA','GBM','HNSC','KIRC','KIRP','LAML','LGG','LIHC',
    #                       'LUAD','LUNG','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')

    #For shortening the list for testing
    #CancerIdentifiers <- c('KICH')

    #These variables are re-defined in the call to the function that computes the Cox model, thus setting these values here does not do anything
    QueryGenes <- c('SLC25A44')
    CoxParameters <- c('Age',QueryGenes) #can also have 'Age' as a CoxParameter

    #Currently plotting is commented out, so setting these values does not do anything
    YRangePlot <- c(-3,3)
    YTicks <- seq(YRangePlot[1],YRangePlot[2],1)

    #create the directory 'output' one level above the current directory for storing the heatmap
    PlotDepositDirectory <- StoreHeatmap()
    DataRepositoryDirectory <- '/Users/Nate/Dropbox/Research/Lehtio_Laboratory/Databases/TCGA_RNA_Seq/' #personal computer
    #DataRepositoryDirectory <- '/Z/users/nate/nate/Data/TCGA_RNA_Seq/' #tamarindo
    #tamarindo working directory command: setwd('/Z/users/nate/nate/functions_and_scripts/heatmap')

    #get a list of all of the directories (one for each tumor-type) that have expression and patient survival data in them
    NumCancers <- length(CancerIdentifiers)
    DataDirectories <- NULL
    for (i in 1:NumCancers)
    {
        if (CancerIdentifiers[i]=='ESCA' || CancerIdentifiers[i]=='STAD') {placeholder <- paste(DataRepositoryDirectory,'TCGA_',CancerIdentifiers[i],'_','exp_HiSeq-2015-02-24/',sep='')}
        else {placeholder <- paste(DataRepositoryDirectory,'TCGA_',CancerIdentifiers[i],'_','exp_HiSeqV2-2015-02-24/',sep='')}
        DataDirectories <- c(DataDirectories,placeholder)
    }

    #name the two files that respectively contain survival data and gene expression data
    PatientInfoFile <- 'clinical_data'
    GeneExpressionFile <- 'genomicMatrix'

    #initialize lists to contain data frames with results for each of the tumor-types
    CoxResultsRisk <- vector('list',NumCancers)

    #find the genes that are common to all cancer-type measurements
    gene_array <- GetCommonGeneList(DataDirectories,GeneExpressionFile)
    gene_array <- sort(gene_array)

    #shorten for testing
    #gene_array <- gene_array[18686:18687]
    #also set the CancerIdentifiers above to a fewer amount if desired

    #iterate through the tumor-types and fit a Cox model to the survival data for each gene
    #for (i in 1:NumCancers)
    for (i in 1:NumCancers)
    {
        print(CancerIdentifiers[i])
        CoxResultsRisk[[i]] <- GetCox(DataDirectories[i],PatientInfoFile,GeneExpressionFile,CoxParameters,CancerIdentifiers[i],gene_array)
    }

    CoxRisk <- do.call("cbind", CoxResultsRisk[1:NumCancers])

#create the directory 'output' one level above the current directory for storing the heatmap
output_directory <- StoreHeatmap()

saveRDS(CoxRisk, file=paste(output_directory,'CoxRisk.rds'))

#Make the heatmap - requires too much memory to be done on a personal computer
hm <- MakeHeatMap(NULL,NULL,NULL,c(-1.5,-0.5,0.5,1.5),NULL,'euclidian','ward.D2',CoxRisk,NULL,NULL,c('SLC25A44','GLS2','GLUL','GOT1','GOT2'),TRUE,FALSE)

return(CoxRisk)
}
