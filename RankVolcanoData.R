#This script creates a descending list of ranked genes based on p-value and fold change of expression between groups

#create a data frame to perform a regression on the columns
#DATA_reg <- DATA_annotated_less_groups #This is the output of the MakeVolcanoPlot function
#genes_to_label <- c('ATIC')

#Regression is performed on a single line, thus the portion of the volcano corresponding to negative regulation (log2(ratio) < 0) is flipped below the y-axis
DATA_reg$lin_nlog10p <- DATA_reg$nlog10p
DATA_reg$lin_nlog10p[DATA_reg$log2ratio<0] <- -DATA_reg$lin_nlog10p[DATA_reg$log2ratio<0] #lin_log10p is the linear regression variable for log10p

#Because the scale of -log10(p) is arbitrary, the scale is normalized to the spread (IQR) in the value
#This is done for the x (log2(ratio)) and y (-log10(p)) axes)
DATA_reg$lin_nlog10p <- DATA_reg$lin_nlog10p/IQR(DATA_reg$lin_nlog10p)
DATA_reg$lin_log2ratio <- DATA_reg$log2ratio/IQR(DATA_reg$log2ratio)

#Perform the linear regression and extract the slope and y-intercept
lin_model <- lm(DATA_reg$lin_nlog10p ~ DATA_reg$lin_log2ratio)
m <- lin_model$coefficients[2]
d <- lin_model$coefficients[1]

#Find the x-coordinates of the intersections points between the lines perpendicular to the best fit line and going through each data point
#The distance of this intersection point along the line of best fit is a measure of how regulated each point (gene) is. The further to the right, the more positively regulated the gene is.
#Since the x-coordinate of the intersection has a positive relationship with the distance along the line, the x-coordinate is satisfactory for ranking regulation
DATA_reg$rank <- (DATA_reg$lin_log2ratio+DATA_reg$lin_nlog10p*m-d*m)/(m^2+1)
DATA_reg <- DATA_reg[order(-DATA_reg$rank),]
ranked_gene_list <- rownames(DATA_reg)
reverse_ranked_gene_list <- rev(ranked_gene_list)

write(ranked_gene_list,paste(volcano_directory,'descending_regulated_gene_list.txt',sep=''))

#Plot the relationship for visualization
#Note the axis must be the same scale for propper visualization of right angles to the fit line
XData <- 'lin_log2ratio' #the name of the column of the data frame to be plotted
YData <- 'lin_nlog10p' #the name of the row of the data frame to be plotted
XAxisLimits <- c(-10*IQR(DATA_reg$lin_log2ratio),10*IQR(DATA_reg$lin_log2ratio))
YAxisLimits <- c(-10*IQR(DATA_reg$lin_nlog10p),10*IQR(DATA_reg$lin_nlog10p))
AxisLabels <- c('Basal/Other Expression','-Log10(p-value)')
PlotDimensions <- c(4,4)

ScatterPlot <-  ggplot(DATA_reg,aes_string(x=XData,y=YData)) +
  #geom_point(color='black') +
  
  #################HERE!!!!!!
  #trying to use the subset which requires the plyr package, need to update other packages and restart r to use it - this is to ensure a selected point is visible on top of the others
  geom_point(color='black') + #specifying one color may remove the legend, May be necessary to make this a string when considering groups
  geom_smooth(method = "lm", se = FALSE, col='red',size=0.1) +  #the fitline?
  annotate('point',DATA_reg[genes_to_label,'lin_log2ratio'],DATA_reg[genes_to_label,'lin_nlog10p'],col='red') + #points to highlight in red on the volcano plot, x-points first, then y-points
  annotate('text',DATA_reg[genes_to_label,'lin_log2ratio'],DATA_reg[genes_to_label,'lin_nlog10p'],col='#0088B2',label=rownames(DATA_reg[genes_to_label,]),size=0.75) + #points to highlight in red on the volcano plot, x-points first, then y-points
  #geom_errorbar(aes_string(ymin=y_error_low_name,ymax=y_error_high_name)) + #adds error bars to y-variable
  #geom_errorbarh(aes_string(xmin=x_error_low_name,xmax=x_error_high_name)) + #adds error bars to x-variable
  #geom_text(label = rownames(DATA_reg),size=0.5,color='red') + #if you want to label the points
  coord_cartesian(xlim=XAxisLimits,ylim=YAxisLimits,expand = 0) + #expand=0 removes extra space before and after first and last tick marks
  #If you need to specify the axis tick marks, uncomment below and comment out coord_cartesian line above
  #scale_x_continuous(limits=XAxisLimits) +
  #scale_y_continuous(limits=YAxisLimits,breaks=c(-4,-3,-2,-1,0,1,2,3,4)) +
  theme(axis.text.y=element_text(color='black',size=10)) +
  theme(axis.ticks.y=element_line(colour='black',size=1)) +
  theme(axis.ticks.x=element_line(colour='black',size=1)) +
  theme(axis.text.x=element_text(color='black',size=10)) +
  theme(axis.title.y=element_text(color='black',vjust=1,size=10)) +
  theme(axis.title.x=element_text(color='black',vjust=0,size=10)) +
  theme(panel.background=element_rect(fill=NA)) +
  theme(axis.line.x=element_line(colour='black',size=1,linetype='solid')) +
  theme(axis.line.y=element_line(colour='black',size=1,linetype='solid')) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +#removes gridlines
  theme(legend.title=element_blank()) +
  theme(legend.key=element_rect(fill=NA)) +
  theme(legend.text = element_text(colour="black", size = 6)) +
  #theme(legend.position=LegendPosition) +
  theme(legend.key.size = unit(0.3, "cm")) + #Space between legend symbols and text, maybe?
  theme(legend.background = element_rect(fill="transparent",linetype = 0,colour = 'transparent')) +
  labs(x=AxisLabels[1]) +
  labs(y=AxisLabels[2]) 
#annotate("text",x=CorrelationLabelPosition[1],y=CorrelationLabelPosition[2],label=paste('PC = ',pearson_coefficient,' [',pearson_95_confidence[1],', ',pearson_95_confidence[2],']','\n','SC = ',spearman_coefficient,' [',spearman_p_value,']',sep=''),col='red',size=2.5,hjust=0,vjust=1)
#annotate("text",x=CorrelationLabelPosition[1],y=CorrelationLabelPosition[2],label=paste('PC = ',pearson_coefficient,' [',pearson_95_confidence[1],', ',pearson_95_confidence[2],']',sep=''),col='red',size=3,hjust=0,vjust=1)

print('work!') #There needs to be something in between ggplot call and ggsave call inside a function. This is some sort of bug in ggplot/ggsave

ggsave(paste(volcano_directory,'volcano_for_regression.pdf',sep=''), width = PlotDimensions[1], height = PlotDimensions[2], dpi = 120)
#rm(DATA_t_annotated) #remove specified global variable so global environment not affected by function. This causes a warning if ColorIndicator is not specified (that's okay)