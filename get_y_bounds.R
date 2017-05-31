get_y_bounds <- function(group_order_original,gene_name,DATA_long)
{
    n_gene <- length(gene_name)
    n_group <- length(group_order_original)

    set_min <- NULL
    set_max <- NULL

    for (i in 1:n_gene)
    {
          for (j in 1:n_group)
          {
              indices <- DATA_long[,'gene']==gene_name[i] & DATA_long[,'group']==group_order_original[j]
              value <- DATA_long[indices,'value']

              q1_m_15iqr <- quantile(value)[2]-1.5*IQR(value)
              min_ge_q1_m_15iqr <- min(value[value >= q1_m_15iqr])

              q3_p_15iqr <- quantile(value)[4]+1.5*IQR(value)
              max_le_q3_p_15iqr <- max(value[value <= q3_p_15iqr])

              if (is.null(set_min) || min_ge_q1_m_15iqr < set_min){set_min <- min_ge_q1_m_15iqr}
              if (is.null(set_max) || max_le_q3_p_15iqr > set_max){set_max <- max_le_q3_p_15iqr}
          }
    }

    y_bounds <- c(set_min,set_max)
    return(y_bounds)
}
