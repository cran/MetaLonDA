#' Visualize Longitudinal Feature
#'
#' It starts by finding by walking up the path until it finds the
#'
#' @param df dataframe has the Count, Group, ID, Time
#' @param text feature name
#' @param group_levels The two level's name
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5 # sample size;
#' n.timepoints = 10 # time point;
#' n.group= 2 # number of group;
#' Group = factor(c(rep(0,n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggretage.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' visualizeFeature(aggretage.df, text = rownames(metalonda_test_data)[1], Group)
#' @export
visualizeFeature = function (df, text, group_levels)
{
  cat("Visualizing Feature = ", text, "\n")
  Count=0;Time=0;ID=0;Group=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  
  p = ggplot(df, aes(Time, Count, colour = Group, group = interaction(Group, ID)))
  p = p + geom_point(size=2) + geom_line(size=1) +  theme_bw() +
    ggtitle(paste("Feature = ", text, sep = "")) +
    scale_colour_manual(values = c("skyblue", "pink"),
                        breaks = c("0", "1"),
                        labels = c(group_levels[1], group_levels[2]))+
    theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
          axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
          legend.text=element_text(size=15, face="plain"), legend.title = element_blank()) +
    theme(legend.position="top") + scale_x_continuous(breaks = waiver()) #min(df$Time):max(df$Time))
  # print(p)
  ggsave(filename=paste("Feature_", text, ".tiff", sep=""), dpi = 300, height = 10, width = 15, units = 'cm')
}



#' Visualize Feature with the fitted Splines
#'
#' Plot the longitudinal features along with the fitted splines
#'
#' @param df dataframe has the Count , Group, ID, Time
#' @param model the fitted model
#' @param text feature name
#' @param method The fitting method (negative binomial, LOWESS)
#' @param group_levels The two levels name
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
visualizeFeatureSpline = function (df, model, method, text, group_levels)
{ 
  cat("Visualizing the Splines of Feature =  ", text, "\n")
  Count=0;Time=0;ID=0;Group=0;lnn=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  ddNULL = model$ddNULL
  dd0 = model$dd0
  dd1 = model$dd1
  ddNULL_U = model$ddNULL_U95
  ddNULL_L = model$ddNULL_L95
  dd0_U = model$dd0_U95
  dd0_L = model$dd0_L95
  dd1_U = model$dd1_U95
  dd1_L = model$dd1_L95
  
  ln = factor(c(rep("longdash", nrow(df)), rep("longdash", nrow(dd0)), rep("longdash", nrow(dd1)), rep("solid", nrow(dd0_U)),
                rep("solid", nrow(dd0_L)), rep("solid", nrow(dd1_U)), rep("solid", nrow(dd1_L))))
  size = c(rep(1, nrow(df)), rep(1, nrow(dd0)), rep(1, nrow(dd1)), rep(1, nrow(dd0_U)),
           rep(1, nrow(dd0_L)), rep(1, nrow(dd1_U)), rep(1, nrow(dd1_L)))
  dm <- rbind(df[,c("Time", "Count", "Group", "ID")], dd0, dd1, dd0_U, dd0_L, dd1_U, dd1_L)
  dm$lnn=ln
  dm$sz= size

  p = ggplot(dm, aes(Time, Count, colour = Group, group = interaction(Group, ID)))
  p = p + theme_bw() + geom_point(size=2) + geom_line(aes(linetype=lnn), size=1) + 
    ggtitle(paste("using ", method,", for Feature = ", text, sep = "")) +
    scale_colour_manual(values = c("skyblue", "pink", "blue", "firebrick",
                                   "blue",  "blue",  "firebrick", "firebrick"), 
                        breaks = c("0", "1", "Fit_0", "Fit_1"),
                        labels = c(group_levels[1], group_levels[2], paste(group_levels[1], "_fit", sep=""), paste(group_levels[2], "_fit", sep="")))+
    theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
          axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
          axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"), 
          legend.text=element_text(size=15, face="plain"), legend.title = element_blank()) +
    theme(legend.position="top") + scale_x_continuous(breaks = waiver()) + guides(linetype=FALSE, size =FALSE)
  
  # print(p)
  ggsave(filename=paste("Feature_", text, " CurveFitting_", method, ".tiff", sep=""), dpi = 300, height = 10, width = 15, units = 'cm')
}



#' Visualize histogram of the Area Ratio
#'
#' Visualize histogram of the Area Ratio emprical distrbition of one time interval
#'
#' @param permuted Permutation of the permuted data
#' @param text Feature name
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
visualizeARHistogram = function(permuted, text){
  cat("Visualizing AR histogram of feature = ", text, "\n")
  tiff(paste("Feature_", text, "_AR_NULL_model.tiff", sep=""), res = 300, height = 20, width = 20, units = 'cm')
  par(mfrow=c(3,3))
  cat("dim of permu = ", dim(permuted), "\n")
  for( i in 1:ncol(permuted)){
    hist(permuted[,i], xlab = "AR Ratio", ylab = "Frequency", 
         breaks = 50, col = "yellow", border = "black", 
         main = paste("AR Null Distribution, interval = ", i, sep=""), xlim = c(0,1))
    lines(density(permuted[,i]), col="blue", lwd = 2)
    lines(density(permuted[,i], adjust=2), lty="dotted", col="darkgreen", lwd=2)
  }
  dev.off()
}



#' Visualize significant time interval
#'
#' It highlights the significant tiem intervals 
#'
#' @param aggretage.df Dataframe has the Count , Group, ID, Time
#' @param model_ss The fitted model
#' @param method Fitting method (negative binomial or LOWESS)
#' @param start Vector of the start points of the time intervals
#' @param end Vector of the end points of the time intervals
#' @param text Feature name
#' @param group_levels Level's name
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
visualizeArea = function(aggretage.df, model_ss, method, start, end, text, group_levels)
{
  cat("Visualizing significant intervals of feature = ", text, "\n")
  Time=0 ## This line is just to pass the CRAN checks for the aes in ggplot2
  sub11 <- list()
  sub10 <- list()
  xx = NULL
  for(i in 1:length(start))
  {
    sub11[[i]] = subset(model_ss$dd1, Time >= start[i] & Time <= end[i])  
    sub10[[i]] = subset(model_ss$dd0, Time >= start[i] & Time <= end[i])
    cmd = sprintf('geom_ribbon(data=sub10[[%d]], aes(ymin=sub11[[%d]]$Count,ymax=Count), colour= "grey3",fill="grey69", alpha="0.6")', i, i)
    if (i != 1)
    {
      xx = paste(xx, cmd, sep = "+")
    } else
    {
      xx = cmd
    }
  }
  
  # ddNULL = model_ss$ddNULL
  dd0 = model_ss$dd0
  dd1 = model_ss$dd1
  
  dm <- rbind(aggretage.df[,c("Time", "Count", "Group", "ID")], dd0, dd1)
  p1 = 'ggplot(dm, aes(Time, Count, colour = Group, group = interaction(Group, ID))) + 
  theme_bw() + geom_point(size = 2) + geom_line(size = 1) + 
  ggtitle(paste("Feature = ", text, sep = "")) +
  scale_colour_manual(values = c("skyblue", "pink", "blue", "firebrick"), 
  breaks = c("0", "1", "Fit_0", "Fit_1"),
  labels = c(group_levels[1], group_levels[2], paste(group_levels[1], "_fit", sep=""), paste(group_levels[2], "_fit", sep=""))) +
  theme(axis.text.x = element_text(colour="black", size=12,angle=0,hjust=0.5,vjust=0.5, face="bold"),
  axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
  axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
  axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"), 
  legend.text=element_text(size=15, face="plain"), legend.title = element_blank()) +
  theme(legend.position="top") + scale_x_continuous(breaks = waiver())' 
  p2 =  xx #"geom_ribbon(data=sub10[[1]], aes(ymin=sub11[[1]]$Count,ymax=Count, group=factor(1), colour = factor(1)), fill=\"green\", alpha=\"0.2\")+geom_ribbon(data=sub10[[2]], aes(ymin=sub11[[2]]$Count,ymax=Count, group=factor(1), colour = factor(1)), fill=\"green\", alpha=\"0.2\")" 
  p3 = paste(p1, p2, sep="+")
  p = eval(parse(text = p3))
  # print(p)
  ggsave(filename=paste("Feature_", text, " SignicantInterval_", method, ".tiff", sep=""), dpi = 300, height = 10, width = 15, units = 'cm')
}

