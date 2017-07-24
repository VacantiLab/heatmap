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
    #CancerIdentifiers <- c('ACC')

    #These variables are re-defined in the call to the function that computes the Cox model, thus setting these values here does not do anything
    QueryGenes <- c('SLC25A44')
    CoxParameters <- c('Age',QueryGenes) #can also have 'Age' as a CoxParameter

    #Currently plotting is commented out, so setting these values does not do anything
    YRangePlot <- c(-3,3)
    YTicks <- seq(YRangePlot[1],YRangePlot[2],1)

    #create the directory 'output' one level above the current directory for storing the heatmap
    PlotDepositDirectory <- StoreHeatmap()
    DataRepositoryDirectory <- '/Users/Nate/Dropbox/Research/Lehtio_Laboratory/Databases/TCGA_RNA_Seq/'

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

    #shorten for testing
    #NumCancers <- 4
    #DataDirectories <- DataDirectories[1:NumCancers]

    #find the genes that are common to all cancer-type measurements
    gene_array <- GetCommonGeneList(DataDirectories,GeneExpressionFile)
    gene_array <- sort(gene_array)

    #shorten for testing
    #gene_array <- gene_array[1:50]

    #iterate through the tumor-types and fit a Cox model to the survival data for each gene
    #for (i in 1:NumCancers)
    for (i in 1:NumCancers)
        {
        print(CancerIdentifiers[i])
        CoxResultsRisk[[i]] <- GetCox(DataDirectories[i],PatientInfoFile,GeneExpressionFile,CoxParameters,CancerIdentifiers[i],gene_array)
        #MakeBoxPlot_survival(BoxPlotDataFrameList[[i]],PlotDepositDirectory,paste(CancerIdentifiers[i],'BoxPlot.pdf',sep='_'),'Factor','Measurement','Risk',YRangePlot,QueryGenes) #the last 3 are x_var, y_var, color_var
        #MakeSurvivalPlot(SurvivalDataFrameList[[i]],PlotDepositDirectory,paste(CancerIdentifiers[i],'SurvivalPlot.pdf',sep='_'))
        #MakeSurvivalCoefficientPlot(CoxResultsDataFrameList[[i]],PlotDepositDirectory,paste(CancerIdentifiers[i],'HazardCoefficientPlot.pdf'),'CoxParameters','CoxParameters') #the last one is the x_var
        }

    CoxRisk <- do.call("cbind", CoxResultsRisk[1:NumCancers])

    ##iterate through the tumor-types and determine if the Cox model is significant
    #CoxModelSignificant <- NULL
    #for (i in 1:NumCancers)
    #{
    #    SignificanceCondition <- CoxResultsDataFrameList[[i]]$CoxLikelihoodRatioWaldP[1]<0.05
    #    CoxModelSignificant <- c(CoxModelSignificant,SignificanceCondition)
    #}

    ##Plotting is not currently functional
    ##Plot the gene expression for each cancer type for each risk group
    ##bind all of the box plot data frames into one data frame
    #AllBoxPlotDataFrame <- do.call('rbind',BoxPlotDataFrameList)
    ##Select only entries for the query gene
    #AllBoxPlotDataFrame <- AllBoxPlotDataFrame[AllBoxPlotDataFrame$Factor == QueryGenes[1],]
    ##select only entries where the separation by risk group and model is significant
    #AllBoxPlotDataFrame <- AllBoxPlotDataFrame[AllBoxPlotDataFrame$RiskLogRankP < 0.05,]
    #AllBoxPlotDataFrame <- AllBoxPlotDataFrame[AllBoxPlotDataFrame$CoxLikelihoodRatioWaldP < 0.05,]
    #AllBoxPlotDataFrame <- AllBoxPlotDataFrame[AllBoxPlotDataFrame$CoxGeneP < 0.05,]
    #AllBoxPlot <- MakeBoxPlot_survival(AllBoxPlotDataFrame,PlotDepositDirectory,'AllBoxPlot.pdf','Identifier','Measurement','Risk',YRangePlot,QueryGenes) #the last 3 are x_var, y_var, color_var

    ##Plot the cox model coefficients
    #AllCoxResultsDataFrame <- do.call('rbind',CoxResultsDataFrameList)
    #AllCoxResultsDataFrame <- AllCoxResultsDataFrame[AllCoxResultsDataFrame$RiskLogRankP < 0.05,] #one measure of significance
    #AllCoxResultsDataFrame <- AllCoxResultsDataFrame[AllCoxResultsDataFrame$CoxLikelihoodRatioWaldP < 0.05,] #another measure of significance
    #AllCoxResultsDataFrame <- AllCoxResultsDataFrame[AllCoxResultsDataFrame$CoxGeneP < 0.05,] #this is not a measure of significance (?), maybe only for the gene - but want age and gene considered
    #AllModelCoefficientPlot <- MakeSurvivalCoefficientPlot(AllCoxResultsDataFrame,PlotDepositDirectory,'AllModelCoefficientPlot.pdf','Identifier','CoxParameters') #the last one is x_var

#survival_ran = TRUE

#create the directory 'output' one level above the current directory for storing the heatmap
output_directory <- StoreHeatmap()

saveRDS(CoxRisk, file=paste(output_directory,'CoxRisk.rds'))
browser()

hm <- MakeHeatMap(NULL,NULL,NULL,c(-1.5,-0.5,0.5,1.5),NULL,'euclidian','ward.D2',CoxRisk,NULL,NULL,c('SLC25A44','GLS2','GLUL','GOT1','GOT2'),TRUE,FALSE)

return(CoxRisk)
}
