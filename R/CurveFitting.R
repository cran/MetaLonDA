#' Fit longitudinal data
#'
#' Fits longitudinal samples from the same group using Negative Binomial or LOWESS
#' 
#' @param df dataframe has the Count, Group, ID, Time
#' @param log log transformation of the data
#' @param method The fitting method (negative binomial, LOWESS)
#' @param points The points at which the prediction should happen
#' @return Returns the fitted model
#' @import gss
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
#' aggretage.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' curveFitting(df = aggretage.df, log = FALSE, method= "nbinomial", points)
#' @export
curveFitting <- function(df, log = FALSE, method = "nbinomial", points){
  if(log == TRUE){
    df$Count =  log2(df$Count+1) 
  }
  
  
  # Seprate the two groups
  group0 = df[df$Group==0, ]
  group1 = df[df$Group==1, ]
  groupNULL = df
  


  ## Fitting 
  if(method == "ss"){
    # null model
    modNULL = ssanova(Count ~ Time, data = groupNULL)
    # full model
    mod0 = ssanova(Count ~ Time, data=group0)
    mod1 = ssanova(Count ~ Time, data=group1)
  }
  else if (method == "lowess"){
    # null model
    modNULL = loess(Count ~ Time, data = groupNULL)
    # full model
    mod0 = loess(Count ~ Time, data = group0)
    mod1 = loess(Count ~ Time, data = group1)
  }
  else if(method == "nbinomial"){
    # null model
    modNULL = gssanova(Count ~ Time, data = groupNULL, family = "nbinomial", skip.iter=TRUE)
    # full model
    mod0 = gssanova(Count ~ Time, data=group0, family = "nbinomial", skip.iter=TRUE)
    mod1 = gssanova(Count ~ Time, data=group1, family = "nbinomial", skip.iter=TRUE)
    mod0_nbinomial_project = project(mod0, c("Time"))
    mod1_nbinomial_project = project(mod1, c("Time"))
    modNULL_nbinomial_project = project(modNULL, c("Time"))
  }
  
  
  ### Calculate goodness of fit F-statistics for the non nbinomial models
  if(method !="nbinomial")
  {
    rssNULL = summary(modNULL)$rss
    rssFULL = summary(mod0)$rss+summary(mod1)$rss
    F_stat = (rssNULL-rssFULL)/rssNULL
  }
  
  
  ## Estimate values at the provided time points
  estNULL = predict(modNULL,data.frame(Time = points), se=TRUE)
  est0 = predict(mod0,data.frame(Time = points), se=TRUE)
  est1 = predict(mod1, data.frame(Time = points), se=TRUE)

  ## prepare dataframe for plotting
  if(method != "nbinomial")
  {
    ### Curve dataframe
    ddNULL = data.frame(Time = points, Count = estNULL$fit, Group = "NULL", ID = "NULL")
    dd0 = data.frame(Time = points, Count = est0$fit, Group = "Fit_0", ID = "Fit_0")
    dd1 = data.frame(Time = points, Count = est1$fit, Group = "Fit_1", ID = "Fit_1")
    
    ### Confidence interval dataframe
    ddNULL_U95 = data.frame(Time = points, Count = (estNULL$fit +1.96*estNULL$se), Group = "NULL_U", ID = "NULL_U")
    ddNULL_L95 = data.frame(Time = points, Count = (estNULL$fit -1.96*estNULL$se), Group = "NULL_L", ID = "NULL_L")
    dd0_U95 = data.frame(Time = points, Count = (est0$fit +1.96*est0$se), Group = "Fit_0_U", ID = "Fit_0_U")
    dd0_L95 = data.frame(Time = points, Count = (est0$fit -1.96*est0$se), Group = "Fit_0_L", ID = "Fit_0_L")
    dd1_U95 = data.frame(Time = points, Count = (est1$fit +1.96*est1$se), Group = "Fit_1_U", ID = "Fit_1_U")
    dd1_L95 = data.frame(Time = points, Count = (est1$fit -1.96*est1$se), Group = "Fit_1_L", ID = "Fit_1_L")
  } else{
    ### Curve dataframe
    ddNULL = data.frame(Time = points, Count = modNULL$nu/exp(estNULL$fit), Group = "NULL", ID = "NULL")
    dd0 = data.frame(Time = points, Count = mod0$nu/exp(est0$fit), Group = "Fit_0", ID = "Fit_0")
    dd1 = data.frame(Time = points, Count = mod1$nu/exp(est1$fit), Group = "Fit_1", ID = "Fit_1")
    
    ### Confidence interval dataframe
    ddNULL_U95 = data.frame(Time = points, Count = modNULL$nu/exp(estNULL$fit +1.96*estNULL$se), Group = "NULL_U", ID = "NULL_U")
    ddNULL_L95 = data.frame(Time = points, Count = modNULL$nu/exp(estNULL$fit -1.96*estNULL$se), Group = "NULL_L", ID = "NULL_L")
    dd0_U95 = data.frame(Time = points, Count = mod0$nu/exp(est0$fit + 1.96*est0$se), Group = "Fit_0_U", ID = "Fit_0_U")
    dd0_L95 = data.frame(Time = points, Count = mod0$nu/exp(est0$fit - 1.96*est0$se), Group = "Fit_0_L", ID = "Fit_0_L")
    dd1_U95 = data.frame(Time = points, Count = mod1$nu/exp(est1$fit + 1.96*est1$se), Group = "Fit_1_U", ID = "Fit_1_U")
    dd1_L95 = data.frame(Time = points, Count = mod1$nu/exp(est1$fit - 1.96*est1$se), Group = "Fit_1_L", ID = "Fit_1_L")
    
  }
  
  ## Return the results
  if(method != "nbinomial")
  {
    output = list(F_stat = F_stat, rssNULL = rssNULL, rssFULL = rssFULL, 
                  ddNULL = ddNULL, dd0 = dd0, dd1 = dd1, modNULL = modNULL, 
                  mod0 = mod0, mod1 = mod1, ddNULL_U95 = ddNULL_U95, ddNULL_L95= ddNULL_L95,
                  dd0_U95 = dd0_U95, dd0_L95= dd0_L95, dd1_U95 = dd1_U95, dd1_L95= dd1_L95)
  }
  else
  {
    output = list(ddNULL = ddNULL, dd0 = dd0, dd1 = dd1, modNULL = modNULL, mod0 = mod0, 
                  mod1 = mod1, mod0_nbinomial_project = mod0_nbinomial_project, 
                  mod1_nbinomial_project = mod1_nbinomial_project, 
                  modNULL_nbinomial_project = modNULL_nbinomial_project,
                  ddNULL_U95 = ddNULL_U95, ddNULL_L95= ddNULL_L95,
                  dd0_U95 = dd0_U95, dd0_L95= dd0_L95, dd1_U95 = dd1_U95, dd1_L95= dd1_L95)
  }
  return(output)
}


#' Calculate Area Ratio of time intervals
#'
#' Fits longitudinal samples from the same group using nbinomial or LOWESS
#' 
#' @param curveFit.df gss data object of the fitted spline
#' @return returns the area ratio for all time intervals
#' @import caTools
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
#' model = curveFitting(df = aggretage.df, log = FALSE, method= "nbinomial", points)
#' intervalArea(model)
#' @export
intervalArea = function(curveFit.df){
  size = length(curveFit.df$ddNULL$Time)
  AR = numeric(size - 1)
  
  ## Calculate the absoulte and the sign of each interval area
  AR_abs = numeric(size - 1)
  AR_sign = numeric(size - 1)
  for(i in 1:(size - 1)){
    area0 = trapz(curveFit.df$dd0$Time[i:(i+1)], curveFit.df$dd0$Count[i:(i+1)])
    area1 = trapz(curveFit.df$dd1$Time[i:(i+1)], curveFit.df$dd1$Count[i:(i+1)])
    areaNULL = trapz(curveFit.df$ddNULL$Time[i:(i+1)], curveFit.df$ddNULL$Count[i:(i+1)])
    AR[i] = (area0 - area1) / max(area0, area1)
    AR_abs[i] = abs(AR[i])
    AR_sign[i] = AR[i]/abs(AR[i])
  }
  
  return(list(AR = AR, AR_abs = AR_abs, AR_sign = AR_sign))
}



#' Find Significant Time Intervals
#'
#' Identify the significant time intervals
#' 
#' @param adjusted_pvalue vector of adjested p-value
#' @param threshold p_value cut off
#' @return returns a list of the start and end points of all significant time intervals
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @examples 
#' p=c(0.04, 0.01, 0.02, 0.04, 0.06, 0.2, 0.06, 0.04)
#' findSigInterval(p, threshold = 0.05)
#' @export
findSigInterval = function(adjusted_pvalue, threshold=0.05)
{
  sig = which(adjusted_pvalue<threshold)
  
  start = numeric()
  end = numeric()
  
  if(length(sig) == 0)
  {
    cat("No significant inteval found \n")
  }
  else if(length(sig) == 1)
  {
    start = sig[1]
    end = sig [1]
  }
  else
  {
    start = sig[1]

    if((sig[2] - sig[1]) != 1)
    {
      end = c(end, sig[1])
    }
    
    for(i in 2:length(sig))
    {
      if(i != length(sig))
      {
        if((sig[i]-sig[i-1]) >1)
        {
          start= c(start, sig[i])
        }
        
        if((sig[i+1] - sig[i]) != 1)
        {
          end = c(end, sig[i])
        }
      }
      else
      {
        if((sig[i]-sig[i-1]) > 1)
        {
          start= c(start, sig[i])
        }
        end= c(end, sig[i])
      }
    }
  }
  
  return(list(start = start, end = end))
}


#' Calculate Area Ratio of time intervals for all permutations
#'
#' Fits longitudinal samples from the same group using negative binomial or LOWESS for all permutations
#' 
#' @param perm list has all the permutated models
#' @return returns a list of all permutation area ratio
#' @references
#' Ahmed Metwally (ametwa2@uic.edu)
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5 # sample size;
#' n.timepoints = 10 # time point;
#' n.perm = 3
#' n.group= 2 # number of group;
#' Group = factor(c(rep(0,n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggretage.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, ID = ID)
#' perm = permutation(aggretage.df, n.perm = 3, method = "nbinomial", points)
#' areaPermutation(perm)
#' @export
areaPermutation = function(perm)
{
  AR_list = list()
  list.len = length(perm)
  for (j in 1:list.len)
  {
    AR_list[[j]] = intervalArea(perm[[j]])
  }
  
  return(AR_list)
}
