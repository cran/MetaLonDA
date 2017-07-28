#' Permute group labels
#'
#' Permutes the group label of the samples in order to construct the empirical distibution 
#'
#' @param permDat dataframe has the Count, Group, ID, Time
#' @param n.perm number of permutations
#' @param method The fitting method (negative binomial, LOWESS)
#' @param points The points at which the prediction should happen
#' @return returns the fitted model for all the permutations
#' @import plyr
#' @import utils
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
#' permutation(aggretage.df, n.perm = 3, method = "nbinomial", points)
#' @export
permutation = function(permDat, n.perm = 10, method = "nbinomial", points){

  n.subjects = length(unique(permDat$ID))
  cat("# of subjects = ", n.subjects, "\n")
  
  ### Solve the problem with this # of combinations
  if(n.subjects < 5){

    PP=list() 
    # for(j in 1:ceiling((n.subjects-1)/2)) ### Why should we divide by 2 then subtract 1?
    # {
      
      comb = combn(1:n.subjects, (n.subjects/2))
      ncomb = dim(comb)[2]  
      PP <- append(PP, llply(1:ncomb, function(i){
        cat("Permutation = ", i, "\n")
        permDat$Group = 0 
        permDat[permDat$ID %in% comb[,i],]$Group = 1 
        perm = curveFitting(df = permDat, log = FALSE, method = method, points)
        assign(paste("Model", i, sep = "_"), perm)
      }, .parallel=FALSE))
    # }
  }
  else
  {
    PP=list() 
    
    PP<- llply(1:n.perm, function(j){
      cat("Permutation = ", j,  "\n")
      for( i in 1: n.subjects){
        perm.uniq.len = 1
        while(perm.uniq.len == 1){
          permDat[which(permDat$ID == i),]$Group = rep(sample(c(0,1),1), sum(permDat$ID == i))# ,  replace = TRUE) #rep(sample(1:2), each = time.point)
          perm.uniq.len = length(unique(permDat$Group))
        }
      }
      perm = curveFitting(df = permDat, log = FALSE, method = method, points)
      assign(paste("Model", j, sep = "_"), perm)
    }, .parallel=FALSE)
  }
  
  return(PP)
}  