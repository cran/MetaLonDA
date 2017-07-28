## ---- fig.show = 'hold', eval = FALSE------------------------------------
#  library(metalonda)
#  data(metalonda_test_data)
#  n.sample = 5 # sample size;
#  n.timepoints = 10 # time point;
#  n.group= 2 # number of group;
#  Group = factor(c(rep(0,n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#  Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#  ID = factor(rep(1:(2*n.sample), each = n.timepoints))
#  points = seq(1, 10, length.out = 10)
#  
#  
#  output_all_nbinomial = metalondaAll(data = metalonda_test_data, Time = Time, Group = Group,
#                                      ID = ID, log = FALSE, fit.method = "nbinomial", n.perm = 10,
#                                      points = points, pvalue_threshold=0.05)
#  
#  output_1_nbinomial = metalonda(Count = metalonda_test_data[1,], Time = Time, Group = Group,
#                                 ID = ID, log = log, fit.method = "nbinomial", n.perm = 10, points = points,
#                                 text=rownames(metalonda_test_data)[1], parall = FALSE, pvalue_threshold=0.05, adjust.method="BH")

