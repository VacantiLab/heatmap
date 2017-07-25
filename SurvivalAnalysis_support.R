CenterScale <- function(DataFrame,Columns,gene)
# This function scales the expression data for survival analysis
# This will be done later by the transform() function
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


GetCox <- function(DataDirectory,PatientInfoFile,GeneExpressionFile,CoxParameters,DataFrameIdentifier,gene_array)
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

    #Extract the required patient information
    PatientInfoRequired <- PatientInfo[,c('X_EVENT','X_OS','days_to_birth')] #days to birth is a negative number, so this is negative age. kept this way because when normalized it is SD below average age
    colnames(PatientInfoRequired) <- c('event','days','Age')

    #Replace the dashes in the Patient Information Parsed Data frame
    SampleIDs <- gsub('-','\\.',PatientInfo$sampleID)
    rownames(PatientInfoRequired) <- SampleIDs

    #Replace the dashes in the genes of the gene_array in order to make valid Cox formula objects with them
    gene_array <- gsub('[[:punct:]]','_',gene_array)
    gene_array <- gsub('^_','X_',gene_array)
    rownames(GeneExpressionInfo) <- gsub('[[:punct:]]','_',rownames(GeneExpressionInfo))
    rownames(GeneExpressionInfo) <- gsub('^_','X_',rownames(GeneExpressionInfo))

    #Determine the number of unique patients
    NumPatients <- length(GeneExpressionInfo[1,])
    NumCoxParameters <- length(CoxParameters)

    #Get a list of the unique patients
    Patients <- colnames(GeneExpressionInfo)

    #Initialize a list of patients whose expression values are missing
    MissingClinicalData <- NULL

    CoxResultsDataFrameList_gene <- rep(list(NULL),length(gene_array))
    names(CoxResultsDataFrameList_gene) <- gene_array

    print_iterator = 1
    for (gene in gene_array)
    {
        print(paste(print_iterator,':',DataFrameIdentifier,':',gene))
        print_iterator = print_iterator+1


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

        #Get the Cox Results for the gene and age
        CoxResultsDataFrame <- GetCoxOutputs(CoxFormula,PatientInfoRequired,gene,CoxParameters,DataFrameIdentifier)

        #Keep the gene results and discard the age results
        CoxResultsDataFrameList_gene[[gene]] <- CoxResultsDataFrame
        CoxResultsDataFrameList_gene[[gene]] <- CoxResultsDataFrameList_gene[[gene]][gene,]

        #remove the gene information from the data frame passed to the CoxModel as it will be repopulated for another gene in the next loop
        PatientInfoRequired[,gene] <- NULL
    }

    #Combine all Dataframes specific for genes by their rows to get a CoxResultsDataFrame for the cancer-type
    CoxResultsDataFrame <- do.call('rbind',CoxResultsDataFrameList_gene)
    CoxResultsDataFrame[,'CoxParameters'] <- NULL
    CoxResultsDataFrame[,'Identifier'] <- NULL

    #Make a data frame to return that only has the risk information
    CoxResultsRisk <- CoxResultsDataFrame[,'Risk',drop=FALSE]
    colnames(CoxResultsRisk) <- DataFrameIdentifier

    return(CoxResultsRisk)
}
