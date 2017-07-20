SurvivalAnalysis <- function()
{
    #load requisite libraries
    library(ggplot2)
    library(survival)
    library(tidyr)
    library(matrixStats)

    CancerIdentifiers <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','COADREAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC',
                           'LUAD','LUNG','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')
    #CancerIdentifiers <- c('ACC')

    QueryGenes <- c('SLC25A44')
    CoxParameters <- c('Age',QueryGenes) #can also have 'Age' as a CoxParameter
    YRangePlot <- c(-3,3)
    YTicks <- seq(YRangePlot[1],YRangePlot[2],1)

    #create the directory 'output' one level above the current directory for storing the heatmap
    PlotDepositDirectory <- StoreHeatmap()
    DataRepositoryDirectory <- '/Users/Nate/Dropbox/Research/Lehtio_Laboratory/Databases/TCGA RNA Seq/'

    NumCancers <- length(CancerIdentifiers)
    DataDirectories <- NULL
    for (i in 1:NumCancers)
    {
        if (CancerIdentifiers[i]=='ESCA' || CancerIdentifiers[i]=='STAD') {placeholder <- paste(DataRepositoryDirectory,'TCGA_',CancerIdentifiers[i],'_','exp_HiSeq-2015-02-24/',sep='')}
        else {placeholder <- paste(DataRepositoryDirectory,'TCGA_',CancerIdentifiers[i],'_','exp_HiSeqV2-2015-02-24/',sep='')}
        DataDirectories <- c(DataDirectories,placeholder)
    }

    PatientInfoFile <- 'clinical_data'
    GeneExpressionFile <- 'genomicMatrix'

    #initialize data frames
    BoxPlotDataFrameList <- vector('list',NumCancers)
    SurvivalDataFrameList <- vector('list',NumCancers)
    CoxResultsDataFrameList <- vector('list',NumCancers)
    KaplanMeierRiskList <- vector('list',NumCancers)

    for (i in 1:NumCancers)
        {
        print(i)
        CoxReturn <- GetCox(DataDirectories[i],PatientInfoFile,GeneExpressionFile,CoxParameters,CancerIdentifiers[i],QueryGenes)
        BoxPlotDataFrameList[[i]] <- CoxReturn[[1]]
        SurvivalDataFrameList[[i]] <- CoxReturn[[2]]
        CoxResultsDataFrameList[[i]] <- CoxReturn[[3]]
        browser()
        #MakeBoxPlot_survival(BoxPlotDataFrameList[[i]],PlotDepositDirectory,paste(CancerIdentifiers[i],'BoxPlot.pdf',sep='_'),'Factor','Measurement','Risk',YRangePlot,QueryGenes) #the last 3 are x_var, y_var, color_var
        #MakeSurvivalPlot(SurvivalDataFrameList[[i]],PlotDepositDirectory,paste(CancerIdentifiers[i],'SurvivalPlot.pdf',sep='_'))
        #MakeSurvivalCoefficientPlot(CoxResultsDataFrameList[[i]],PlotDepositDirectory,paste(CancerIdentifiers[i],'HazardCoefficientPlot.pdf'),'CoxParameters','CoxParameters') #the last one is the x_var
        }

    CoxModelSignificant <- NULL
    for (i in 1:NumCancers)
    {
        SignificanceCondition <- CoxResultsDataFrameList[[i]]$CoxLikelihoodRatioWaldP[1]<0.05
        CoxModelSignificant <- c(CoxModelSignificant,SignificanceCondition)
    }

    #Plot the gene expression for each cancer type for each risk group
    #bind all of the box plot data frames into one data frame
    AllBoxPlotDataFrame <- do.call('rbind',BoxPlotDataFrameList)
    #Select only entries for the query gene
    AllBoxPlotDataFrame <- AllBoxPlotDataFrame[AllBoxPlotDataFrame$Factor == QueryGenes[1],]
    #select only entries where the separation by risk group and model is significant
    #AllBoxPlotDataFrame <- AllBoxPlotDataFrame[AllBoxPlotDataFrame$RiskLogRankP < 0.05,]
    #AllBoxPlotDataFrame <- AllBoxPlotDataFrame[AllBoxPlotDataFrame$CoxLikelihoodRatioWaldP < 0.05,]
    #AllBoxPlotDataFrame <- AllBoxPlotDataFrame[AllBoxPlotDataFrame$CoxGeneP < 0.05,]
    AllBoxPlot <- MakeBoxPlot_survival(AllBoxPlotDataFrame,PlotDepositDirectory,'AllBoxPlot.pdf','Identifier','Measurement','Risk',YRangePlot,QueryGenes) #the last 3 are x_var, y_var, color_var

    #Plot the cox model coefficients
    AllCoxResultsDataFrame <- do.call('rbind',CoxResultsDataFrameList)
    #AllCoxResultsDataFrame <- AllCoxResultsDataFrame[AllCoxResultsDataFrame$RiskLogRankP < 0.05,] #one measure of significance
    #AllCoxResultsDataFrame <- AllCoxResultsDataFrame[AllCoxResultsDataFrame$CoxLikelihoodRatioWaldP < 0.05,] #another measure of significance
    #AllCoxResultsDataFrame <- AllCoxResultsDataFrame[AllCoxResultsDataFrame$CoxGeneP < 0.05,] #this is not a measure of significance (?), maybe only for the gene - but want age and gene considered
    AllModelCoefficientPlot <- MakeSurvivalCoefficientPlot(AllCoxResultsDataFrame,PlotDepositDirectory,'AllModelCoefficientPlot.pdf','Identifier','CoxParameters') #the last one is x_var

survival_ran = TRUE
return(survival_ran)
}
