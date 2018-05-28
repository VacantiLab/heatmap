transform_data <- function(DATA,transformation)
{
  transformed = FALSE

  if (is.null(transformation)){return(DATA)}

  #Log transform the data if specified to do so
  if (transformation == 'log2')
  {
    DATA <- log2(DATA)
    #DATA[DATA < -3] <- -3
    transformed = TRUE
  }

    #Log10 transform the data if specified to do so
    if (transformation == 'log10')
    {
      DATA <- log10(DATA)
      #DATA[DATA < -3] <- -3
      transformed = TRUE
    }

  #exp2() transform the data if specified to do so
  if (transformation == 'exp2')
  {
    DATA <- 2^(DATA)
    transformed = TRUE
  }

  #center rows on median and normalize by IQR
  if (transformation == 'median_center_iqr_norm')
  {
    column_names <- colnames(DATA) #record the column names after the unecessary column is removed
    Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
    DATA <- data.frame(lapply(Transposed_DATA, median_center_iqr_norm))
    DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
    colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
    transformed = TRUE
  }

  #raise to power 2, median norm, then log2 transform
  if (transformation == 'exp2_mednorm_log2')
  {
    column_names <- colnames(DATA) #record the column names after the unecessary column is removed
    Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
    DATA <- data.frame(lapply(Transposed_DATA, exp2_mednorm_log2))
    DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
    colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
    transformed = TRUE
  }

  #Normalize rows to row median and log2 transform
  if (transformation == 'median_norm_log2_transform')
  {
    column_names <- colnames(DATA) #record the column names after the unecessary column is removed
    Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
    DATA <- data.frame(lapply(Transposed_DATA, median_norm_log2_transform))
    DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
    colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
    transformed = TRUE
  }

  #normalize by iqr, log2 trasform, then center on resulting median
  if (transformation == 'iqr_norm_log2_median_center')
  {
    column_names <- colnames(DATA) #record the column names after the unecessary column is removed
    Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
    DATA <- data.frame(lapply(Transposed_DATA, iqr_norm_log2_median_center))
    DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
    colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
    transformed = TRUE
  }

  #raise to power 2, normalize by iqr, log2 transform, then center on resulting median
  if (transformation == 'exp2_iqrnorm_log2_mediancenter')
  {
    column_names <- colnames(DATA) #record the column names after the unecessary column is removed
    Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
    DATA <- data.frame(lapply(Transposed_DATA, exp2_iqrnorm_log2_mediancenter))
    DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
    colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
    transformed = TRUE
  }

  if (transformed==FALSE && !is.null(transformation)){stop('custom message: You have specified a transformation that does not exist.')}

  return(DATA)
}


###############################################################################


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

#Function to normalize by iqr, log2 trasform, then center on resulting median
iqr_norm_log2_median_center <- function(vector)
{
  #normalized <- vector/median(vector)
  normalized <- vector/IQR(vector)
  transformed <- log2(normalized)
  transformed <- transformed - median(transformed)
  return(transformed)
}


#Function to raise to power 2, median norm, then log2 transform
exp2_mednorm_log2 <- function(vector)
{
  exp2_raised <- vector^2
  normalized <- exp2_raised/median(exp2_raised)
  transformed <- log2(normalized)
  return(transformed)
}


#Function to raise to power 2, iqr norm, log2 transform, then median center
exp2_iqrnorm_log2_mediancenter <- function(vector)
{
  exp2_raised <- vector^2
  normalized <- exp2_raised/IQR(exp2_raised)
  transformed <- log2(normalized)
  transformed <- transformed - median(transformed)
  return(transformed)
}
