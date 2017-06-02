get_y_bounds <- function(group_order,gene_name,DATA_long)
# This function examines each box of the box and whisker plot and returns the y-bounds of the graph
# These bounds correspond to the lowest point of all of the whiskers and the highest point of all of the whiskers
{
    n_gene <- length(gene_name)
    n_group <- length(group_order)
    if (is.null(group_order)){n_group=1}

    set_min <- NULL
    set_max <- NULL

    for (i in 1:n_gene)
    {
          for (j in 1:n_group)
          {
              #For each gene-group combination, get arrays of values
              if (!is.null(group_order)){indices <- DATA_long[,'gene']==gene_name[i] & DATA_long[,'group']==group_order[j]}
              if (is.null(group_order)){indices <- DATA_long[,'gene']==gene_name[i]}
              value <- DATA_long[indices,'value']

              #Find value corresponding to the bottom of the whisker - the smallest value >= 25th percentile - 1.5*IQR
              q1_m_15iqr <- quantile(value)[2]-1.5*IQR(value)
              min_ge_q1_m_15iqr <- min(value[value >= q1_m_15iqr])

              #Find value corresponding to the top of the whisker - the largest value <= 75th percentile + 1.5*IQR
              q3_p_15iqr <- quantile(value)[4]+1.5*IQR(value)
              max_le_q3_p_15iqr <- max(value[value <= q3_p_15iqr])

              #Record the lowest and highest of these values respectively for all gene-group combinations
              if (is.null(set_min) || min_ge_q1_m_15iqr < set_min){set_min <- min_ge_q1_m_15iqr}
              if (is.null(set_max) || max_le_q3_p_15iqr > set_max){set_max <- max_le_q3_p_15iqr}
          }
    }

    #Return the y-bounds
    y_bounds <- c(set_min,set_max)
    return(y_bounds)
}
