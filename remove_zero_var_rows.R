remove_zero_var_rows <- function(DATA)
{
  var_vec <- apply(DATA,1,var)
  non_zero_var_indices <- which(var_vec!=0)
  DATA <- DATA[non_zero_var_indices,]
}