RankVolcanoData <- function(volcano_df,output_directory)
{
#This script creates a descending list of ranked genes based on p-value and fold change of expression between groups

#create a data frame to perform a regression on the columns
#DATA_reg <- DATA_annotated_less_groups #This is the output of the MakeVolcanoPlot function
#genes_to_label <- c('ATIC')

#Regression is performed on a single line, thus the portion of the volcano corresponding to negative regulation (log2(ratio) < 0) is flipped below the y-axis
volcano_df$reg_nlog10_p <- volcano_df$nlog10_p
volcano_df$reg_nlog10_p[volcano_df$log2_ratio<0] <- -volcano_df$reg_nlog10_p[volcano_df$log2_ratio<0] #lin_log10p is the linear regression variable for log10p

#Because the scale of -log10(p) is arbitrary, the scale is normalized to the spread (IQR) in the value
#This is done for the x (log2(ratio)) and y (-log10(p)) axes)
volcano_df$reg_nlog10_p <- volcano_df$reg_nlog10_p/IQR(volcano_df$reg_nlog10_p)
volcano_df$reg_log2_ratio <- volcano_df$log2_ratio/IQR(volcano_df$log2_ratio)

#Perform the linear regression and extract the slope and y-intercept
lin_model <- lm(volcano_df$reg_nlog10_p ~ volcano_df$reg_log2_ratio)
m <- lin_model$coefficients[2]
d <- lin_model$coefficients[1]

#Find the x-coordinates of the intersections points between the lines perpendicular to the best fit line and going through each data point
#The distance of this intersection point along the line of best fit is a measure of how regulated each point (gene) is. The further to the right, the more positively regulated the gene is.
#Since the x-coordinate of the intersection has a positive relationship with the distance along the line, the x-coordinate is satisfactory for ranking regulation
volcano_df$rank <- (volcano_df$reg_log2_ratio+volcano_df$reg_nlog10_p*m-d*m)/(m^2+1)
volcano_df <- volcano_df[order(-volcano_df$rank),]
ranked_gene_list <- rownames(volcano_df)

write(ranked_gene_list,paste(output_directory,'descending_regulated_gene_list.txt',sep=''))

return(volcano_df)
}
