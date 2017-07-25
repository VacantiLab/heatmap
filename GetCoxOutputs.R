GetCoxOutputs <- function(CoxFormula,PatientInfoRequired,gene,CoxParameters,DataFrameIdentifier)
# Sometimes the coxph will throw up an error because the data doesn't quite work
#     In these cases the risk for that gene and that cancer-type is returned as 0 and the function continues
#     One instance of thses cases is in KICH for gene TRIM60
{
    #If there is no error thrown by coxph()
    attempt <- try(CoxModel <- coxph(CoxFormula, data = PatientInfoRequired))
    if (!("try-error" %in% class(attempt)))
    {
        CoxCoefficients <- CoxModel$coefficients
        HazardModelCoefficients <- exp(CoxCoefficients) #This are the outputs of the model! Being one SD from the mean for the given measurement multiplies the hazard function by this coefficient
        CoxCoefficientBounds <- confint(CoxModel)
        #     2.5%   97.5%
        #Age
        #Gene

        WaldCoxFitP <- 1-pchisq(CoxModel$wald.test,length(CoxModel$coefficients)) #the p for the whole model
        CoxGeneP <- summary(CoxModel)$coefficients[gene,5] #The p for the gene coefficient

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
        CoxResultsDataFrame$Risk <- NULL
        if (is.na(CoxResultsDataFrame[gene,'HazardModelCoefficientLowerBounds'])){CoxResultsDataFrame[gene,'Risk'] <- 0}
        if (!is.na(CoxResultsDataFrame[gene,'HazardModelCoefficientLowerBounds']))
        {
            if (CoxResultsDataFrame[gene,'HazardModelCoefficientLowerBounds'] < 1 && CoxResultsDataFrame[gene,'HazardModelCoefficientUpperBounds'] > 1){CoxResultsDataFrame[gene,'Risk'] <- 0}
            if (CoxResultsDataFrame[gene,'HazardModelCoefficientLowerBounds'] > 1){CoxResultsDataFrame[gene,'Risk'] <- 1}
            if (CoxResultsDataFrame[gene,'HazardModelCoefficientUpperBounds'] < 1){CoxResultsDataFrame[gene,'Risk'] <- -1}
        }
    }

    #If there is an error thrown by coxph() return the risk as 0
    #The other fields do not matter now because they are not considered in the analysis
    if ("try-error" %in% class(attempt))
    {
        CoxResultsDataFrame <- data.frame(matrix(nrow=2,ncol=5))
        rownames(CoxResultsDataFrame) <- c(gene,'Age')
        colnames(CoxResultsDataFrame) <- c('CoxParameters','Identifier','CoxLikelihoodRatioWaldP','CoxGeneP','Risk')
        CoxResultsDataFrame$CoxParameters <- c(0,0)
        CoxResultsDataFrame$Identifier <- c(0,0)
        #CoxResultsDataFrame$RiskLogRankP <- rep(RiskLogRankP,length(CoxResultsDataFrame[,1]))
        CoxResultsDataFrame$CoxLikelihoodRatioWaldP <- c(0,0)
        CoxResultsDataFrame$CoxGeneP <- c(0,0)
        CoxResultsDataFrame$Risk <- c(0,0)
    }

    return(CoxResultsDataFrame)
}
