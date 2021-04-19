transform_data <- function(DATA,transformation,select_rows_after_transform)
{
    transformed = FALSE
    if (is.null(transformation))
    {
        transformation <- 'none'
        transformed <- TRUE
    }

    #Log transform the data if specified to do so
    if (transformation == 'log2')
    {
        #First replace 0 values with half the row minimum
        nCols <- ncol(DATA)
        DATA[1:nCols] <- lapply(DATA[1:nCols], function(x) replace(x, x == 0, min(x[x>0]/2)))
        #log2 transform
        DATA <- log2(DATA)
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
    if (transformation == 'row_median_center_iqr_norm')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, median_center_iqr_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    #center rows on mean and normalize by sample SD
    if (transformation == 'row_mean_center_sd_norm')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, mean_center_sd_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    #log2, center rows on median, and normalize by IQR
    if (transformation == 'log2_row_median_center_iqr_norm')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        DATA <- log2(DATA)
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, median_center_iqr_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    #center rows on median and normalize by IQR
    if (transformation == 'row_median_center_range_norm')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, median_center_range_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'row_median_center_iqr_norm_col_iqr_norm')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, median_center_iqr_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function

        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        DATA <- data.frame(lapply(DATA, iqr_norm))
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names

        transformed = TRUE
    }

    #raise to power 2, median norm, then log2 transform
    if (transformation == 'exp2_row_mednorm_log2')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, exp2_mednorm_log2))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    #Normalize rows to row median and log2 transform
    if (transformation == 'row_median_norm_log2')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, median_norm_log2_transform))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    #Normalize rows to row average and log2 transform
    if (transformation == 'row_average_norm_log2')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, average_norm_log2_transform))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    #normalize by iqr, log2 trasform, then center on resulting median
    if (transformation == 'row_iqr_norm_log2_median_center')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, iqr_norm_log2_median_center))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    #raise to power 2, normalize by iqr, log2 transform, then center on resulting median
    if (transformation == 'row_exp2_iqrnorm_log2_mediancenter')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, exp2_iqrnorm_log2_mediancenter))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'row_mednorm_col_mednorm_log2')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, median_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        DATA <- data.frame(lapply(DATA, median_norm_log2_transform)) #median norm the columns and then log2 transform everything
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'row_mednorm_col_mednorm')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        Transposed_DATA <- data.frame(t(DATA)) #normalize to row median, transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, median_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        DATA <- data.frame(lapply(DATA, median_norm)) #median norm the columns
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'row_meannorm_col_mednorm_log2')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, mean_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        DATA <- data.frame(lapply(DATA, median_norm_log2_transform)) #median norm the columns and then log2 transform everything
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'add_halfnon0min_row_meannorm_col_mednorm_log2')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        DATA <- DATA + 0.5*min(DATA[DATA!=0])
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, mean_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        DATA <- data.frame(lapply(DATA, median_norm_log2_transform)) #median norm the columns and then log2 transform everything
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'add_halfnon0min_row_meannorm_col_mednorm')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        DATA <- DATA + 0.5*min(DATA[DATA!=0])
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, mean_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        DATA <- data.frame(lapply(DATA, median_norm)) #median norm the columns and then log2 transform everything
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'add_halfnon0min_row_meannorm_log2')
    {
      column_names <- colnames(DATA) #record the column names after the unecessary column is removed
      row_names <- rownames(DATA)
      DATA <- DATA + 0.5*min(DATA[DATA!=0])
      Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
      DATA <- data.frame(lapply(Transposed_DATA, mean_norm))
      DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
      DATA <- log2(DATA)
      colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
      rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
      transformed = TRUE
    }

    if (transformation == 'add_halfnon0min_row_meannorm')
    {
      column_names <- colnames(DATA) #record the column names after the unecessary column is removed
      row_names <- rownames(DATA)
      DATA <- DATA + 0.5*min(DATA[DATA!=0])
      Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
      DATA <- data.frame(lapply(Transposed_DATA, mean_norm))
      DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
      colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
      rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
      transformed = TRUE
    }

    if (transformation == 'zrow')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, zscore))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'log2_zrow')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        DATA <- log2(DATA)
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, zscore))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformation == 'row_meannorm_col_mednorm_log2_zrow')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, mean_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        DATA <- data.frame(lapply(DATA, median_norm_log2_transform)) #median norm the columns and then log2 transform everything
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, zscore)) #zscore transform rows
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }
    
    if (transformation == 'row_meannorm_col_mednorm_log2_zrow_zcol')
    {
      column_names <- colnames(DATA) #record the column names after the unecessary column is removed
      row_names <- rownames(DATA)
      Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
      DATA <- data.frame(lapply(Transposed_DATA, mean_norm))
      DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
      DATA <- data.frame(lapply(DATA, median_norm_log2_transform)) #median norm the columns and then log2 transform everything
      Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
      DATA <- data.frame(lapply(Transposed_DATA, zscore)) #zscore transform rows
      DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
      DATA <- data.frame(lapply(DATA, zscore)) #zscore transform columns
      colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
      rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
      transformed = TRUE
    }

    if (transformation == 'add_halfnon0min_row_meannorm_col_mednorm_log2_zrow')
    {
        column_names <- colnames(DATA) #record the column names after the unecessary column is removed
        row_names <- rownames(DATA)
        DATA <- DATA + 0.5*min(DATA[DATA!=0])
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, mean_norm))
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        DATA <- data.frame(lapply(DATA, median_norm_log2_transform)) #median norm the columns and then log2 transform everything
        Transposed_DATA <- data.frame(t(DATA)) #transpose because can only scale columns
        DATA <- data.frame(lapply(Transposed_DATA, zscore)) #zscore transform rows
        DATA <- data.frame(t(DATA)) #transpose back resulting in scaled rows
        colnames(DATA) <- column_names #give the column names back because they are lost when converted to a matrix by t() function
        rownames(DATA) <- row_names #give the row names back because they are lost when converted to a matrix by t() function
        transformed = TRUE
    }

    if (transformed==FALSE && !is.null(transformation)){stop('custom message: You have specified a transformation that does not exist.')}

    if (!is.null(select_rows_after_transform))
    {
        select_rows_after_transform_in_data_set <- select_rows_after_transform %in% rownames(DATA)
        select_rows_after_transform <- select_rows_after_transform[select_rows_after_transform_in_data_set]
        DATA <- DATA[select_rows_after_transform,]
    }

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

mean_center_sd_norm <- function(vector)
{
  centered <- vector - mean(vector)
  scaled <- centered/sd(vector)
  return(scaled)
}

#Function to center on median and scale by IQR
median_center_range_norm <- function(vector)
{
  centered <- vector - median(vector)
  scaled <- centered/(max(vector)-min(vector))
  return(scaled)
}

#Function to center on median and scale by IQR
iqr_norm <- function(vector)
{
  centered <- vector
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

#Function to normalize by median and log2 transform
average_norm_log2_transform <- function(vector)
{
  normalized <- vector/mean(vector)
  transformed <- log2(normalized)
  return(transformed)
}

#Function to normalize by median and log2 transform
median_norm <- function(vector)
{
  normalized <- vector/median(vector)
  transformed <- normalized
  return(transformed)
}

#Function to normalize by median and log2 transform
mean_norm <- function(vector)
{
  normalized <- vector/mean(vector)
  transformed <- normalized
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

#Function to transform to z scores
zscore <- function(vector)
{
  transformed <- (vector - mean(vector))/sd(vector)
  return(transformed)
}
