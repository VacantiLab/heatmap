transform_data <- function(DATA,transformation)
{
  transformed = FALSE

  #Log transform the data if specified to do so
  if (transformation == 'log2')
  {
    DATA <- log2(DATA)
    DATA[DATA < -3] <- -3
    transformed = TRUE
  }

  #Center and scale the rows if indicated
  if (transformation == 'median_center_iqr_norm')
  {
    column_names <- colnames(DATA) #record the column names after the unecessary column is removed
    Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
    DATA <- data.frame(lapply(Transposed_DATA, median_center_iqr_norm))
    DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
    colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
    transformed = TRUE
  }

  #Normalize rows to row median and log2 transform if indicated
  if (transformation == 'median_norm_log2_transform')
  {
    column_names <- colnames(DATA) #record the column names after the unecessary column is removed
    Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
    DATA <- data.frame(lapply(Transposed_DATA, median_norm_log2_transform))
    DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
    colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
    transformed = TRUE
  }

  if (transformed==FALSE && !is.null(transformation)){stop('custom message: You have specified a transformation that does not exist.')}

  return(DATA)
}

#Function to center on median and scale by IQR
median_center_iqr_norm <- function(vector)
{
  centered <- vector - median(vector)
  scaled <- centered/IQR(vector)
  return(scaled)
}

#Function to normalize by median and log2 transform
median_norm_log2_transform <- function(vector)
{
  normalized <- vector/median(vector)
  transformed <- log2(normalized)
  return(transformed)
}
