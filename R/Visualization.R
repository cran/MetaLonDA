#' Visualize Longitudinal Feature
#'
#' Visualize Longitudinal Feature
#'
#' @param df dataframe has the Count, Group, ID, Time
#' @param text feature name
#' @param group.levels The two level's name
#' @param unit time interval unit
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures (default: "Normalized Count")
#' @param prefix prefix to be used to create directory for the analysis results
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' data(metalonda_test_data)
#' pfx = tempfile()
#' dir.create(file.path(pfx))
#' n.sample = 5
#' n.timepoints = 10
#' n.group = 2
#' Group = factor(c(rep(0, n.sample*n.timepoints), rep(1, n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggregate.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' visualizeFeature(df = aggregate.df, text = rownames(metalonda_test_data)[1], 
#' group.levels = Group, prefix = pfx)
#' @export
visualizeFeature = function (df, text, group.levels, unit = "days", ylabel = "Normalized Count", 
                             col = c("blue", "firebrick"), prefix = "Test")
{
  cat("Visualizing Feature = ", text, "\n")
  Count=0; Time=0; ID=0; Group=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  
  p = ggplot(df, aes(Time, Count, colour = Group, group = interaction(Group, ID)))
  p = p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
    ggtitle(paste("Feature = ", text, sep = "")) + labs(y = ylabel, x = sprintf("Time (%s)", unit)) +
    scale_colour_manual(values = col, breaks = c("0", "1"),
                        labels = c(group.levels[1], group.levels[2])) +
    theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
          axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
          legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
          plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="top") + scale_x_continuous(breaks = waiver())
  
  #print("Prefix = ", prefix)
  #ggsave(filename=paste("Feature_", text, ".jpg", sep=""), dpi = 1200, height = 10, width = 15, units = 'cm')
  ggsave(filename=paste(prefix, "/", "Feature_", text, ".jpg", sep=""), dpi = 1200, height = 10, width = 15, units = 'cm')
}



#' Visualize the feature trajectory with the fitted Splines
#'
#' Plot the longitudinal features along with the fitted splines
#'
#' @param df dataframe has the Count , Group, ID, Time
#' @param model the fitted model
#' @param method The fitting method (negative binomial, LOWESS)
#' @param group.levels The two level's name
#' @param text feature name
#' @param unit time unit used in the Time vector (hours, days, weeks, months, etc.)
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures (default: "Normalized Count")
#' @param prefix prefix to be used to create directory for the analysis results
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
visualizeFeatureSpline = function (df, model, method, text, group.levels, unit = "days", ylabel = "Normalized Count", 
                                   col = c("blue", "firebrick"), prefix = "Test")
{ 
  cat("Visualizing Splines of Feature = ", text, "\n")
    
  Count=0;Time=0;ID=0;Group=0;lnn=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  dd.null = model$dd.null
  dd.0 = model$dd.0
  dd.1 = model$dd.1
  
  ln = factor(c(rep("longdash", nrow(df)), rep("longdash", nrow(dd.0)), rep("longdash", nrow(dd.1))))
  size = c(rep(1, nrow(df)), rep(1, nrow(dd.0)), rep(1, nrow(dd.1)))
  dm = rbind(df[,c("Time", "Count", "Group", "ID")], dd.0, dd.1)
  dm$lnn=ln
  dm$sz= size
  
  p = ggplot(dm, aes(Time, Count, colour = Group, group = interaction(Group, ID)))
  p = p + theme_bw() + geom_point(size=1, alpha=0.5) + geom_line(aes(linetype=lnn), size=1, alpha=0.5) + 
    ggtitle(paste("Feature = ", text, sep = "")) + labs(y = ylabel, x = sprintf("Time (%s)", unit)) +
    scale_colour_manual(values = c(col, col), 
                        breaks = c("0", "1", "fit.0", "fit.1"),
                        labels = c(group.levels[1], group.levels[2], paste(group.levels[1], ".fit", sep=""), paste(group.levels[2], ".fit", sep="")))+
    theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
          axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"), 
          legend.text=element_text(size=15, face="plain"), legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="top") + scale_x_continuous(breaks = waiver()) + guides(linetype=FALSE, size =FALSE)
  
  ggsave(filename=paste(prefix, "/", "Feature_", text, "_CurveFitting_", method, ".jpg", sep=""), dpi = 1200, height = 10, width = 15, units = 'cm')
}




#' Visualize significant time interval
#'
#' Visualize significant time interval
#'
#' @param aggregate.df Dataframe has the Count, Group, ID, Time
#' @param model.ss The fitted model
#' @param method Fitting method (negative binomial or LOWESS)
#' @param start Vector of the start points of the time intervals
#' @param end Vector of the end points of the time intervals
#' @param text Feature name
#' @param group.levels Level's name
#' @param unit time unit used in the Time vector (hours, days, weeks, months, etc.)
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures (default: "Normalized Count")
#' @param prefix prefix to be used to create directory for the analysis results
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
visualizeArea = function(aggregate.df, model.ss, method, start, end, text, group.levels, unit = "days", 
                         ylabel = "Normalized Count", col = c("blue", "firebrick"), prefix = "Test")
{
  cat("Visualizing Significant Intervals of Feature = ", text, "\n")
  Time = 0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  sub.11 = list()
  sub.10 = list()
  xx = NULL
  for(i in 1:length(start))
  {
    sub.11[[i]] = subset(model.ss$dd.1, Time >= start[i] & Time <= end[i])  
    sub.10[[i]] = subset(model.ss$dd.0, Time >= start[i] & Time <= end[i])
    cmd = sprintf('geom_ribbon(data=sub.10[[%d]], aes(ymin = sub.11[[%d]]$Count, ymax = Count), colour= "grey3", fill="grey69", 
                  alpha = 0.6)', i, i)
    if (i != 1)
    {
      xx = paste(xx, cmd, sep = "+")
    } else
    {
      xx = cmd
    }
  }
  
  # ddNULL = model_ss$ddNULL
  dd.0 = model.ss$dd.0
  dd.1 = model.ss$dd.1
  
  dm = rbind(dd.0, dd.1)
  p1 = 'ggplot(dm, aes(Time, Count, colour = Group, group = interaction(Group, ID))) + 
  theme_bw() + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.5) + 
  ggtitle(paste("Feature = ", text, sep = "")) + labs(y = ylabel, x = sprintf("Time (%s)", unit)) +
  scale_colour_manual(values = col, 
  breaks = c("fit.0", "fit.1"),
  labels = c(paste(group.levels[1], ".fit", sep = ""), paste(group.levels[2], ".fit", sep = ""))) +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(colour = "black", size = 12, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
  axis.title.x = element_text(colour = "black", size = 15, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
  axis.title.y = element_text(colour = "black", size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"), 
  legend.text = element_text(size = 15, face="plain"), legend.title = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "top") + scale_x_continuous(breaks = waiver())' 
  p2 = xx  
  p3 = paste(p1, p2, sep="+")
  p = eval(parse(text = p3))
  ggsave(filename=paste(prefix, "/", "Feature_", text, "_SignificantInterval_", method, ".jpg", sep=""), dpi = 1200, height = 10, width = 15, units = 'cm')
}



#' Visualize all significant time intervals for all tested features
#'
#' Visualize all significant time intervals for all tested features
#'
#' @param interval.details Dataframe has infomation about significant interval (feature name, start, end, dominant, p-value)
#' @param prefix prefix for the output figure
#' @param unit time unit used in the Time vector (hours, days, weeks, months, etc.)
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
visualizeTimeIntervals = function(interval.details, prefix = "Test", unit = "days", 
                                  col = c("blue", "firebrick"))
{
  feature=0;dominant=0;ID=0;Group=0;lnn=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  interval.details$dominant = as.factor(interval.details$dominant)
  interval.details$pvalue = as.numeric((interval.details$pvalue))
  interval.details = interval.details[order(interval.details$feature), ]
  
  
  ggplot(interval.details, aes(ymin = start , ymax = end, x = feature, xend = feature)) + 
    geom_linerange(aes(color = dominant), size = 1) + 
    coord_flip() +  scale_colour_manual(values = col) +
    labs(x = "Feature", y = sprintf("Time (%s)", unit), colour="Dominant") + 
     theme(axis.text.x = element_text(colour = "black", size = 10, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
           axis.text.y = element_text(colour = "black", size = 8, angle = 0, vjust = 0.5, face = "bold"),
           axis.title.x = element_text(colour = "black", size = 15, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
           axis.title.y = element_text(colour = "black", size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),
           legend.text = element_text(size = 15, face = "plain")) + 
    theme(panel.grid.minor =   element_blank(),
          panel.grid.major.y = element_line(colour = "white", size = 6),
          panel.grid.major.x = element_line(colour = "white",size = 0.75)) +
    theme(legend.position="top", panel.border = element_rect(colour = "black", fill = NA, size = 2))
  ggsave(filename = paste(prefix, "/", prefix, "_MetaLonDA_TimeIntervals.jpg", sep=""), dpi = 1200, height = 30, width = 20, units = 'cm')
}




#' Visualize Area Ratio (AR) empirical distribution
#'
#' Visualize Area Ratio (AR) empirical distribution for each time interval
#'
#' @param permuted Permutation of the permuted data
#' @param text Feature name
#' @param method fitting method
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @param prefix prefix to be used to create directory for the analysis results
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
visualizeARHistogram = function(permuted, text, method, prefix = "Test"){
  cat("Visualizing AR Distribution for Feature = ", text, "\n")
  n = ncol(permuted)
  r = ceiling(sqrt(n))
  c = ceiling(sqrt(n))
  xx = paste(prefix, "/", "Feature_", text, "_AR_distribution_", method, ".jpg", sep = "")
  jpeg(filename = xx, res = 1200, height = r*5, width = c*5, units = 'cm')
  
  par(mfrow=c(r,c))
  for( i in 1:ncol(permuted)){
    hist(permuted[,i], xlab = "AR Ratio", ylab = "Frequency", 
         breaks = 10, col = "yellow", border = "red", 
         main = paste("Interval # ", i, sep=""), xlim = c(0,1))
  }
  dev.off()
}




#' Visualize log2 fold-change and significance of each interval as volcano plot
#'
#' Visualize log2 fold-change and significance of each interval as volcano plot
#'
#' @param df Dataframe has a detailed summary about feature's significant intervals
#' @param text Feature name
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @param prefix prefix to be used to create directory for the analysis results
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
visualizeVolcanoPlot = function(df, text, prefix = "Test"){
  adjusted.pvalue_pseudo=0; Significance=0; log2FoldChange=0 ## This line is just to pass the CRAN checks
  cat("Visualizing Volcano Plot of Feature = ", text, "\n")
  # Highlight features that have an absolute log2 fold change > 1 and p-value < 0.05
  df$adjusted.pvalue_pseudo = df$adjusted.pvalue
  df$adjusted.pvalue_pseudo[which(df$adjusted.pvalue == 0)] = 0.00001
  df$Significance <- "NS"
  df$Significance[(abs(df$log2FoldChange) > 1)] <- "FC"
  df$Significance[(df$adjusted.pvalue_pseudo<0.05)] <- "FDR"
  df$Significance[(df$adjusted.pvalue_pseudo<0.05) & (abs(df$log2FoldChange)>1)] <- "FC_FDR"
  table(df$Significance)
  df$Significance <- factor(df$Significance, levels=c("NS", "FC", "FDR", "FC_FDR"))
  
  
  ggplot(data=df, aes(x=log2FoldChange, y=-log10(adjusted.pvalue_pseudo), colour=Significance)) +
    geom_point(alpha=0.4, size=1.75) + theme_bw() +
    ggtitle(paste("Feature = ", text, sep = "")) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
          axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
          legend.text=element_text(size=10, face="plain"), legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5), legend.position="bottom") +
    
    scale_color_manual(values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), 
                       labels=c(NS="NS", FC=paste("Log2FC>|", 1, "|", sep=""), FDR=paste("FDR Q<", 0.05, sep=""), 
                                FC_FDR=paste("FDR Q<", 0.05, " & Log2FC>|", 1, "|", sep=""))) +
    
    #Tidy the text labels for a subset of genes
    geom_text(data=subset(df, adjusted.pvalue_pseudo<0.05 & abs(log2FoldChange)>=1),
              aes(label=rownames(subset(df, adjusted.pvalue_pseudo<0.05 & abs(log2FoldChange)>= 1))),
              size=2.25,
              #segment.color="black", #This and the next parameter spread out the labels and join them to their points by a line
              #segment.size=0.01,
              check_overlap=TRUE,
              vjust=1.0)+
    geom_vline(xintercept=c(-1,1), linetype="dotted") + geom_hline(yintercept=-log10(0.05), linetype="dotted") +
    ggsave(filename = paste(prefix, "/", "Feature_", text, "_VolcanoPlot.jpg", sep=""), dpi = 1200, height = 12, width = 12, units = 'cm')
}