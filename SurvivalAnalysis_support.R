CenterScale <- function(DataFrame,Columns,gene)
{
  FactorMean <- NULL
  FactorSD <- NULL
  NumColumns <- length(Columns)
  NumRows <- length(DataFrame[,1])
  for (i in 1:NumColumns)
  {
    FactorMean <- c(FactorMean,mean(as.numeric(DataFrame[,Columns[i]]),na.rm=TRUE)) #make this mean!!!!!
    FactorSD_append <- sd(as.numeric(DataFrame[,Columns[i]]),na.rm=TRUE)
    if (FactorSD_append==0){FactorSD_append = 1} #cannot have NaN resulting from division by zero below in Cox Model
    FactorSD <- c(FactorSD,FactorSD_append)
  }
  FactorMeanMatrix <- t(replicate(NumRows, FactorMean))
  FactorSDMatrix <- t(replicate(NumRows, FactorSD))
  ReturnDataFrame <- (DataFrame[Columns]-FactorMeanMatrix)/FactorSDMatrix
  return(ReturnDataFrame)
}
##################################################################################################
GetCox <- function(DataDirectory,PatientInfoFile,GeneExpressionFile,CoxParameters,DataFrameIdentifier,QueryGenes)
{
    # QueryGenes is now defined below to loop through all genes
    # CoxParameters is also defined in the loop because it is looping through Cox models for each gene paired with Age
    # This function is being made to return the following for a cancer type:
    #        HazardCoefficient    UpperBound    LowerBound
    # Gene1
    # Gene2
    #
    # The hazard coefficient is for a model on a single gene and the age of the patient
    # This is being done so that the hazard coefficient profile of genes across all of the cancers can be compared to each other
    #    This will allow the hazard coefficient profiles of each gene to be gathered

    PatientInfoFilePath <- paste(DataDirectory,PatientInfoFile,sep='')
    PatientInfo <- read.table(file=PatientInfoFilePath,head=TRUE,sep='\t')

    GeneExpressionFilePath <- paste(DataDirectory,GeneExpressionFile,sep='')
    GeneExpressionInfo <- read.table(file=GeneExpressionFilePath,head=TRUE,sep='\t')
    rownames(GeneExpressionInfo) <- GeneExpressionInfo$sample
    GeneExpressionInfo$sample <- NULL #Remove the column containing the gene gene names as they are no longer needed
    gene_array <- rownames(GeneExpressionInfo)[1:1000]

    #Extract the required patient information
    PatientInfoRequired <- PatientInfo[,c('X_EVENT','X_OS','days_to_birth')] #days to birth is a negative number, so this is negative age. kept this way because when normalized it is SD below average age
    colnames(PatientInfoRequired) <- c('event','days','Age')

    #Replace the dashes in the Patient Information Parsed Data frame
    SampleIDs <- gsub('-','\\.',PatientInfo$sampleID)
    rownames(PatientInfoRequired) <- SampleIDs

    #Determine the number of unique patients
    NumPatients <- length(GeneExpressionInfo[1,])
    NumQueryGenes <- length(QueryGenes)
    NumCoxParameters <- length(CoxParameters)

    #Get a list of the unique patients
    Patients <- colnames(GeneExpressionInfo)

    #Initialize a list of patients whose expression values are missing
    MissingClinicalData <- NULL

    CoxResultsDataFrameList_gene <- rep(list(NULL),length(gene_array))
    names(CoxResultsDataFrameList_gene) <- gene_array

    for (gene in gene_array)
    {
        print(gene)

        #Initialize the Column for the Gene Expression
        NoEntryIdentifier = 1000000
        PatientInfoRequired[,gene] <- NoEntryIdentifier

        CoxParameters <- c('Age',gene)


        #Fill in the patient expression values
        for (i in 1:NumPatients)
        {
            if (Patients[i] %in% rownames(PatientInfoRequired)){PatientInfoRequired[Patients[i],gene] <- GeneExpressionInfo[gene,Patients[i]]}
            else{MissingClinicalData <- c(MissingClinicalData,Patients[i])}
        }

        #Get a list of the patient sample IDs, note there can be multiple samples per patient
        PatientSampleIDs <- rownames(PatientInfoRequired)
        NumPatientSampleIDs <- length(PatientSampleIDs)

        #Only keep the patient sample IDs corresponding to those samples with expression data
        PatientSampleIDsKeep <- NULL
        for (i in 1:NumPatientSampleIDs)
        {
            if (PatientInfoRequired[PatientSampleIDs[i],gene] != NoEntryIdentifier){PatientSampleIDsKeep <- c(PatientSampleIDsKeep,PatientSampleIDs[i])}
        }
        PatientInfoRequired <- PatientInfoRequired[PatientSampleIDsKeep,]

        #Create the Cox Model
        PatientInfoRequired$SurvObj <- with(PatientInfoRequired,Surv(days,event==1))
        # km.as.one <- survfit(SurvObj ~ 1, data = PatientInformationRequired, conf.type = "log-log")
        # plot(km.as.one)

        #Center and scale the parameter variables
        PatientInfoRequired[CoxParameters] <- CenterScale(PatientInfoRequired,CoxParameters,gene)
        #PatientInfoRequired$Age <- -PatientInfoRequired$Age #Being one means one SD below the average!!! This should be protective!

        #create the Cox model object
        CoxParametersSumString <- paste(CoxParameters,collapse=' + ')
        CoxFormulaString <- paste('SurvObj ~ ',CoxParametersSumString)
        CoxFormula <- as.formula(CoxFormulaString)
        CoxModel <- coxph(CoxFormula, data = PatientInfoRequired)
        WaldCoxFitP <- 1-pchisq(CoxModel$wald.test,length(CoxModel$coefficients)) #the p for the whole model
        CoxGeneP <- summary(CoxModel)$coefficients[gene,5] #The p for the gene coefficient

        ##PI is the risk index
        #CoxCoefficientMatrix <- t(replicate(NumPatients, CoxModel$coefficients))
        #if (length(CoxParameters) > 1){PatientInfoRequired$PI <- rowSums(PatientInfoRequired[,CoxParameters]*CoxCoefficientMatrix)}
        #if (length(CoxParameters) == 1){PatientInfoRequired$PI <- PatientInfoRequired[,CoxParameters]*CoxModel$coefficients[1]}
        #PatientInfoRequired <- PatientInfoRequired[order(PatientInfoRequired[,'PI']),]

        ##Divide into risk groups
        #NumLowRisk <- ceiling(NumPatients/2)
        #NumHighRisk <- NumPatients-NumLowRisk
        #PatientInfoRequired$Risk <- 'Empty'
        #PatientInfoRequired$Risk[1:NumLowRisk] <- 'Low Risk'
        #PatientInfoRequired$Risk[(NumLowRisk+1):NumPatients] <- 'High Risk'
        #PatientInfoRequired$CoxLikelihoodRatioWaldP <- rep(WaldCoxFitP,NumPatients)
        #PatientInfoRequired$CoxGeneP <- rep(CoxGeneP,NumPatients)

        ##Create the Kaplan Meier objects for plotting survival over time for risk groups
        #KaplanMeierRisk <- survfit(SurvObj ~ Risk, data = PatientInfoRequired, conf.type = "log-log")
        #KaplanMeierDiff <- survdiff(SurvObj ~ Risk, data = PatientInfoRequired)
        #RiskLogRankP <- 1-pchisq(KaplanMeierDiff$chisq,length(KaplanMeierDiff$n)-1) #the p for the separate risk groups
        #PatientInfoRequired$RiskLogRankP <- rep(RiskLogRankP,NumPatients)

        #BoxPlotDataFrame <- PatientInfoRequired[c(CoxParameters,'Risk','CoxLikelihoodRatioWaldP','RiskLogRankP','CoxGeneP')]
        #BoxPlotDataFrame <- BoxPlotDataFrame[complete.cases(BoxPlotDataFrame),]
        #BoxPlotDataFrame <- gather_(BoxPlotDataFrame,'Factor','Measurement',CoxParameters)
        #BoxPlotDataFrame$Identifier <- rep(DataFrameIdentifier,length(BoxPlotDataFrame[,1]))
        #BoxPlotDataFrame$Risk <- factor(BoxPlotDataFrame$Risk, levels=c("Low Risk", "High Risk")) #sets the order for the boxplot

        ##Create a merged data frame for low and high risk groups for plotting
        #LowRiskIndicator <- rep('High Risk',KaplanMeierRisk$strata[1]) #I don't know what dermines the order, but based on the numbers in each group this is the correct assignment
        #HighRiskIndicator <- rep('Low Risk',KaplanMeierRisk$strata[2])
        #MergedTime <- KaplanMeierRisk$time
        #MergedSurvival <- KaplanMeierRisk$surv
        #MergedRiskIndicator <- c(LowRiskIndicator,HighRiskIndicator)
        #SurvivalDataFrame <- data.frame(MergedTime,MergedSurvival,MergedRiskIndicator)
        #SurvivalDataFrame$MergedRiskIndicator <- factor(SurvivalDataFrame$MergedRiskIndicator, levels=c("Low Risk", "High Risk")) #sets the order for the boxplot

        #Get the model outputs
        #A cox model is the hazard function multiplied by exp(CoxCoef*risk_factor), or H*exp(Cox*R)
        #The Cox coefficient is the coefficient to the hazard function when R=1 or exp(Cox*1)
        #In this model R is defined as standard deviation above the average for expression and the negative of age
        CoxCoefficients <- CoxModel$coefficients
        HazardModelCoefficients <- exp(CoxCoefficients) #This are the outputs of the model! Being one SD from the mean for the given measurement multiplies the hazard function by this coefficient
        CoxCoefficientBounds <- confint(CoxModel)
        #     2.5%   97.5%
        #Age
        #Gene
        CoxCoefficientLowerBounds <- CoxCoefficientBounds[,1]
        CoxCoefficientUpperBounds <- CoxCoefficientBounds[,2]
        HazardModelCoefficientLowerBounds <- exp(CoxCoefficientLowerBounds)
        HazardModelCoefficientUpperBounds <- exp(CoxCoefficientUpperBounds)

        CoxResultsDataFrame <- data.frame(CoxParameters,HazardModelCoefficients,HazardModelCoefficientLowerBounds,HazardModelCoefficientUpperBounds)
        CoxResultsDataFrame$CoxParameters <- factor(CoxResultsDataFrame$CoxParameters)
        CoxResultsDataFrame$Identifier <- rep(DataFrameIdentifier,length(CoxResultsDataFrame[,1]))
        #CoxResultsDataFrame$RiskLogRankP <- rep(RiskLogRankP,length(CoxResultsDataFrame[,1]))
        CoxResultsDataFrame$CoxLikelihoodRatioWaldP <- rep(WaldCoxFitP,length(CoxResultsDataFrame[,1]))
        CoxResultsDataFrame$CoxGeneP <- rep(CoxGeneP,length(CoxResultsDataFrame[,1]))

        #CoxResultsDataFrame['Age',] <- NULL
        CoxResultsDataFrameList_gene[[gene]] <- CoxResultsDataFrame
        CoxResultsDataFrameList_gene[[gene]] <- CoxResultsDataFrameList_gene[[gene]][gene,]

        PatientInfoRequired[,gene] <- NULL
    }

    CoxResultsDataFrame <- do.call('rbind',CoxResultsDataFrameList_gene)
    CoxResultsDataFrame[,'CoxParameters'] <- NULL
    CoxResultsDataFrame[,'Identifier'] <- NULL

    browser()

    CoxReturn <- list(BoxPlotDataFrame,SurvivalDataFrame,CoxResultsDataFrame)
    return(CoxReturn)
}
################################################################################################
MakeBoxPlot_survival <- function(BoxPlotDataFrame,FileDirectory,FileName,x_var,y_var,color_var,YRangePlot,QueryGenes)
{
  FillColors <- c('#87CEFA','#6495ED')
  names(FillColors) <- c('High Risk','Low Risk')
  TextSize = 12
  if (x_var=='Identifier') {YLabel = paste(QueryGenes[1],' Expression\n','(Standard Deviations From Mean)',sep='')}
  else {YLabel = 'Standard Deviations From Mean'}

  sp <- ggplot(BoxPlotDataFrame,aes_string(x=x_var, y=y_var)) +
        #stat_boxplot(geom ='errorbar') +
        geom_boxplot(aes_string(fill=color_var),outlier.colour='black',outlier.size=0,width=0.7,position=position_dodge(width=0.6)) +
        #scale_y_continuous(breaks=YTicks, limits=YRangePlot, expand=c(0,0)) +
        #coord_cartesian(ylim=c(-3,4)) +
        #stat_summary(position=position_dodge(width=0.5)) +
        #coord_cartesian(ylim=YRangePlot) +
        theme(axis.text.y=element_text(color='black',size=TextSize)) +
        theme(axis.ticks.y=element_line(colour='black',size=1)) +
        theme(axis.ticks.x=element_line(colour='black',size=1)) +
        theme(axis.text.x=element_text(color='black',size=TextSize,angle=90,hjust=1,vjust=0.5)) +
        theme(axis.title.y=element_text(color='black',vjust=1,size=TextSize)) +
        theme(axis.title.x=element_text(color='black',vjust=0,size=TextSize)) +
        theme(panel.background=element_rect(fill=NA)) +
        theme(axis.line=element_line(colour='black',size=1,linetype='solid')) +
        theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + #removes gridlines
        theme(legend.title=element_blank()) +
        theme(legend.key=element_rect(fill=NA)) + #No background color in legend
        theme(legend.text = element_text(colour="black", size=TextSize)) +
        theme(legend.position=c(0.60,0.90)) +
        scale_fill_manual(values=FillColors) +
        #theme(legend.key.size = unit(0.2, "cm")) +
        labs(x = '') +
        labs(y =YLabel)
        #aes_string() allows the factors to be specified by strings and ensures they are evaluated within the correct environment (aes() causes all sorts of trouble)

  ggsave(paste(FileDirectory,FileName,sep=''), width = 12, height = 4, dpi = 120)
  return()
}
########################################################################################################################################
MakeSurvivalPlot <- function(SurvivalDataFrame,FileDirectory,FileName)
{
  FillColors <- c('#87CEFA','#6495ED')
  names(FillColors) <- c('High Risk','Low Risk')
  TextSize = 12

  ggplot(SurvivalDataFrame,aes_string(x='MergedTime',y='MergedSurvival',group='MergedRiskIndicator')) +
    geom_line(aes_string(colour='MergedRiskIndicator')) +
    geom_point(aes_string(colour='MergedRiskIndicator'),size=3) +
    scale_colour_manual(values=FillColors) +
    scale_y_continuous(limits=c(0.2,1),labels=scales::percent,expand=c(0,0)) + #expand=c(0,0) ensures there is no space between data and axis
    scale_x_continuous(limits=c(0,4000),breaks=c(0,1000,2000,3000,4000),expand=c(0,0)) +
    theme(axis.text.y=element_text(color='black',size=TextSize)) +
    theme(axis.ticks.y=element_line(colour='black',size=1)) +
    theme(axis.ticks.x=element_line(colour='black',size=1)) +
    theme(axis.text.x=element_text(color='black',size=TextSize)) +
    theme(axis.title.y=element_text(color='black',vjust=1,size=TextSize)) +
    theme(axis.title.x=element_text(color='black',vjust=0,size=TextSize)) +
    theme(panel.background=element_rect(fill=NA)) +
    theme(axis.line=element_line(colour='black',size=1,linetype='solid')) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + #removes gridlines
    theme(legend.title=element_blank()) +
    theme(legend.key=element_rect(fill=NA)) + #No background color in legend
    theme(legend.text = element_text(colour="black", size=TextSize)) +
    theme(legend.position=c(0.8,0.85)) +
    labs(x = 'Days') +
    labs(y ='Percent Surviving')

  ggsave(paste(FileDirectory,FileName,sep=''), width = 4, height = 4, dpi = 120)
}
###############################################################################################################
MakeSurvivalCoefficientPlot <- function(CoxResultsDataFrame,FileDirectory,FileName,x_var,color_var)
{
  FillColors <- c('#87CEFA','#6495ED')
  TextSize = 12

  ggplot(CoxResultsDataFrame, aes_string(x=x_var,y='HazardModelCoefficients',fill=color_var)) +
    geom_bar(position=position_dodge(),stat='identity',color='black') +
    geom_errorbar(aes_string(ymax='HazardModelCoefficientUpperBounds',ymin='HazardModelCoefficientLowerBounds'),
                  position=position_dodge(0.9),width=0.2,col='black') +
    geom_hline(yintercept=1,linetype=2) +
    #scale_y_continuous(expand=c(0,0),limits=c(0,4.0)) +
    coord_cartesian(ylim=c(0,7)) +
    theme(axis.text.y=element_text(color='black',size=TextSize)) +
    theme(axis.ticks.y=element_line(colour='black',size=1)) +
    theme(axis.ticks.x=element_line(colour='black',size=1)) +
    theme(axis.text.x=element_text(color='black',size=TextSize,angle=90,hjust=1,vjust=0.5)) +
    theme(axis.title.y=element_text(color='black',vjust=1,size=TextSize)) +
    theme(axis.title.x=element_text(color='black',vjust=0,size=TextSize)) +
    theme(panel.background=element_rect(fill=NA)) +
    theme(axis.line=element_line(colour='black',size=1,linetype='solid')) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + #removes gridlines
    theme(legend.title=element_blank()) +
    theme(legend.key=element_rect(fill=NA)) + #No background color in legend
    theme(legend.text = element_text(colour="black", size=TextSize)) +
    theme(legend.position=c(0.40,0.80)) +
    theme(legend.direction='horizontal') +
    scale_fill_manual(values=FillColors) +
    guides(fill = guide_legend(override.aes = list(colour = NULL))) + #removes diagonal line from legend entries
    theme(legend.key = element_rect(colour = 'black')) + #adds border to legened entries
    #theme(legend.key.size = unit(0.2, "cm")) +
    labs(x = '') +
    labs(y ='Hazard Model Coefficient')

  ggsave(paste(FileDirectory,FileName,sep=''), width = 12.0, height = 4, dpi = 120)

}
