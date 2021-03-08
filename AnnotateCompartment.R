AnnotateCompartment <- function(GeneRowDF,CompartmentFile,AnnotationColumnName1='Compartment',AnnotationColumnName2='Neighborhood')
{
  GeneAnnotationDF <- read_txt_to_df(CompartmentFile)
  
  GenesInData <- rownames(GeneRowDF)
  GenesWithCompartment <- rownames(GeneAnnotationDF)
  
  GenesInDataWithCompartmentIndices <- GenesInData %in% GenesWithCompartment
  GenesInDataWithCompartment <- GenesInData[GenesInDataWithCompartmentIndices]
  
  for (gene in GenesInData)
  {
    if (gene %in% GenesInDataWithCompartment)
    {
      GeneRowDF[gene,'Compartment'] <- as.character(GeneAnnotationDF[gene,AnnotationColumnName1])
      GeneRowDF[gene,'Neighborhood'] <- as.character(GeneAnnotationDF[gene,AnnotationColumnName2])
    }
    
    if (!(gene %in% GenesInDataWithCompartment))
    {
      GeneRowDF[gene,'Compartment'] <- 'NotDetected'
      GeneRowDF[gene,'Neighborhood'] <- 'NotDetected'
    }
  }
  
  return(GeneRowDF)

    
}