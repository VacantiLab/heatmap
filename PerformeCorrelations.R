PerformCorrelations <- function(DATA,GeneToCorrelate)
{
  nGenes <- dim(DATA)[1]
  GeneToCorrelateExpression <- DATA[GeneToCorrelate,]
  Genes <- rownames(DATA)

  # Initialize a data frame to hold the pearson correlations, confidence intervals, and p-values
  Correlations <- data.frame(matrix(nrow=nGenes, ncol=5))
  rownames(Correlations) <- rownames(DATA)
  colnames(Correlations) <- c('pearson','U95','L95','p','log10p')
  
  # Calculate correlations, confidence intervals, and p-values and store them in the data frame
  i = 0
  for (Gene in Genes)
  {
    x1 <- as.numeric(DATA[Gene,])
    x2 <- as.numeric(GeneToCorrelateExpression)
    if (sum(x1)==0){Correlations <- Correlations[rownames(Correlations)!=Gene,]}
    if (sum(x1)!=0)
      {
      CorTest <- cor.test(x1, x2, alternative='two.sided', method ="pearson")
      Correlations[Gene,] <- c(CorTest$estimate,CorTest$conf.int[1],CorTest$conf.int[2],CorTest$p.value,-log(CorTest$p.value,10))
      }
    i = i+1
    if (i%%1000==0){print(paste('correlating gene number',toString(i)))}
    #if (i%%1000==0){break}
  }
  # Remove GeneToCorrelateExpression because it, of course, has a perfect correlation with itself
  GeneToCorrelate_row_index <- rownames(Correlations) == GeneToCorrelate
  Correlations <- Correlations[!GeneToCorrelate_row_index,]
  
  
  return(Correlations)


}