#' Metagenomic Longitudinal Differential Abundant Analysis for one feature
#'
#' Find significant time interval of the tested feature
#' 
#' @param Count Count of feature for all groups for all time points for all samples.
#' @param Time Time label of all samples.
#' @param Group Group label of all samples.
#' @param ID individual ID label for samples.
#' @param n.perm number of permutations.
#' @param log log transformation of the data.
#' @param fit.method The fitting method (NB, LOWESS).
#' @param points The points at which the prediction should happen.
#' @param text Feature's name.
#' @param parall Logic to indicate whether to use multicore.
#' @param pvalue_threshold p_value threshold cutoff.
#' @param adjust.method multiple testing correction methods.
#' @return Returns a list of the significant time intervals for the tested feature.
#' @import parallel
#' @import doParallel
#' @import stats
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
#' output_1_nbinomial = metalonda(Count = metalonda_test_data[1,], Time = Time, Group = Group,
#' ID = ID, log = log, fit.method =  "nbinomial", n.perm = 10, points = points,
#' text=rownames(metalonda_test_data)[1], parall = FALSE, pvalue_threshold=0.05, adjust.method="BH")
#' @export
metalonda = function(Count, Time, Group, ID, 
                      n.perm = 10, log = FALSE, fit.method = "nbinomial", 
                      points, text = 0, parall = FALSE, pvalue_threshold=0.05, 
                      adjust.method = "BH")
{
  cat("Start MetaLonDA \n")
  
  group_levels = sort(unique(Group))
  if(length(group_levels) > 2){
    stop("You have more than two phenotypes.")
  }
  gr1 = group_levels[1]
  gr2 = group_levels[2]
  Group[which(Group == gr1)] = 0
  Group[which(Group == gr2)] = 1
  
  
  aggretage.df = data.frame(Count = Count, Time = Time, Group = Group, ID = ID)

  # cat("Visualize feature's trajectories \n")
  visualizeFeature(aggretage.df, text, group_levels)


  
  group0 = aggretage.df[aggretage.df$Group == 0, ]
  group1 = aggretage.df[aggretage.df$Group == 1, ]
  points_min = max(sort(group0$Time)[1], sort(group1$Time)[1])
  points_max = min(sort(group0$Time)[length(group0$Time)], sort(group1$Time)[length(group1$Time)])
  points = points[which(points>=points_min & points<=points_max)]
  

  cat("Start Curve Fitting \n") 
  if (fit.method == "ss")
  {
    cat("Fitting: Gaussian SS \n") 
    model= curveFitting(df = aggretage.df, log = FALSE, method= "ss", points)
  }
  
  else if (fit.method == "lowess")
  {
    cat("Fitting: LOWESS \n")
    model = curveFitting(df = aggretage.df, log = FALSE, method= "lowess", points)
  }

  else if (fit.method == "nbinomial")
  {
    cat("Fitting: NB SS \n")
    model = tryCatch(curveFitting(df = aggretage.df, log = FALSE, method= "nbinomial", points),
             error = function(e) {print(paste("ERROR in gss")); return("ERROR")})
    # if(model=="ERROR")
    #   return("ERROR")
  }
  
  # cat("Visualize feature's trajectories spline \n")
  # visualizeFeatureSpline(aggretage.df, model, fit.method, text, group_levels)
 
   
  ## Calculate area under the curve for each time interval
  cat("Calculate Are of each Unit Interval \n")
  area = intervalArea(model)
  
  
  
  ### Run in Parallel
  if(parall == TRUE) {
    max.cores = detectCores()
    desired.cores = max.cores - 1		
    cl = makeCluster(desired.cores)
    registerDoParallel(cl)
  } 

  
  ## Permutation 
  cat("Start Permutation \n")
  perm  = permutation(aggretage.df, n.perm, fit.method, points)
  

  ### Area P-value per unit interval
  area_perm = areaPermutation(perm)

  if(parall == TRUE) {
    stopCluster(cl)
  } 
  
  
  a1 = do.call(rbind, area_perm)
  a2 = do.call(rbind, a1[,2])

  ## Histogram for one area interval
  ## visualizeARHistogram(a2, text)
  
  ## Calculate area p-value 
  pvalue_area = sapply(1:(length(points)-1), function(i){
    sum(a2[,i] >= area$AR_abs[i])/length(a2[,i])
  } )

  # cat("pvalue_area= ")
  # print(pvalue_area)
  # cat("\n")

  ## Visualize sigificant area
	cat("P_value adjustment method = ", adjust.method, "\n")
  adjusted_pvalue = p.adjust(pvalue_area, method = adjust.method)
  interval = findSigInterval(adjusted_pvalue, threshold = pvalue_threshold)

  st = points[interval$start]
  en = points[interval$end + 1]
  
  if(length(st)>0)
  {
    visualizeArea(aggretage.df, model, fit.method, st, en, text, group_levels)
  }
  
  cat("\n\n")
  
  output_details = list(feature = text, significant_interval=cbind(start = st, end = en), 
       intervals_pvalue=pvalue_area, adjusted_pvalue = adjusted_pvalue, areaSign = area$AR_sign)
  output_summary = data.frame(feature=rep(output_details$feature, nrow(output_details$significant_interval)), 
                 output_details$significant_interval, 
                 dominant = output_details$areaSign[output_details$significant_interval[,1]])
  
  return(list(detailed = output_details, summary = output_summary))
}


#' Metagenomic Longitudinal Differential Abundant Analysis for all features
#'
#' Find significant features and their time interval
#' 
#' @param data Count matrix of all features
#' @param Time Time label of all samples
#' @param Group Group label of all samples
#' @param ID individual ID label for samples
#' @param n.perm number of permutations
#' @param log log transformation of the data
#' @param fit.method The fitting method (NB, LOWESS)
#' @param points The points at which the prediction should happen
#' @param parall logic to indicate whether to use multicore
#' @param pvalue_threshold p-value threshold cutoff
#' @param adjust.method Multiple testing correction methods
#' @return Returns a list of the significant features a long with their significant time intervals
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @export
metalondaAll = function(data, Time, Group, ID, n.perm = 10, log = FALSE, 
                         fit.method = "nbinomial", points, parall = FALSE, 
                         pvalue_threshold =0.05, adjust.method = "BH")
{
  n.features = nrow(data)
  detailed = list()
  summary = list()
  for (i in 1:n.features)
  {
    cat ("Feature  = ", rownames(data)[i], "\n")
    out = metalonda(Count = data[i,], Time = Time, Group = Group, ID = ID, log = log, 
                          fit.method = fit.method, n.perm = n.perm, points = points, 
                          text=rownames(data)[i], parall = parall, pvalue_threshold, adjust.method)
    
    detailed[[i]] = out$detailed
    summary[[i]] = out$summary
  }
  
  summary_tmp = do.call(rbind, summary)
  return(list(output_detail = detailed, output_summary = summary_tmp))
}
