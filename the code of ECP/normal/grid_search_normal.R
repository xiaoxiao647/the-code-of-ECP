##########################################grid_search_normal#################################################

#Generate simulated external data
n0 <- 50  # Sample size of the external control group
u0 <- 1   # True mean value for the external control group
sig0 <- 1   # True variance value for the external control group

##Set random seeds to ensure reproducible results
set.seed(8172)  
x0 <- rnorm(n0, u0, sig0)
u_0 <- mean(x0)
sd_0 <- sd(x0)


## control and treatment
nc <- 25 
nt <- 50

uc <- u_0
ut <- u_0

sigc <- sd_0
sigt <- sd_0


ntrial <- 1000     #Set the number of simulations

gamma <- 1         #Set the value of gamma


#############################################q1#######################
q1 <- seq(0.01, 0.99, by = 0.01) ##
q2 <-  0.01

cl<- makeCluster(23)      
registerDoParallel(cl)   

##Choose q1,q2 to maximize Utility
timestart <- Sys.time()

pp_q1 <- foreach(i = 1:length(q1),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma"))  %dopar% {
  
  #calculate parameter a and b
  para <- decide_para(x0, n0, nc, gamma, q1[i], q2, small = 0.00000000001, large = 100000000000, R = 50000)
  a <- para[1] 
  b <- para[2]
  ##probability that treatment is superior to control
  
  pp <- ECP_normal_model(a, b, n0, x0, nc, uc=u_0, sigc=sd_0, nt, ut=u_0, sigt=sd_0, ntrial)
  return(pp)
}
stopCluster(cl)
timeend <- Sys.time()
print(timeend - timestart)


#Calibrate Cutoff
cutoff_values <- seq(0.001, 0.999, by = 0.001)           #Set the range of cutoff values of probability that treatment is superior to control

selected_cutoff_list <- numeric(length(q1))
#function of calculating probability of rejecting null hypothesis
for (j in 1:length(q1)) {
  rej_null_p <- function(pp, cutoff_values, ntrial){
    prob.rej <- numeric(length(cutoff_values))
    for (k in 1:length(cutoff_values)) {
      cutoff <- cutoff_values[k]
      rej_null <- 0
      
      for(l in 1:ntrial){
        if(pp[l] > cutoff){
          rej_null <- rej_null + 1
        }
      }
      prob.rej[k] <- rej_null/ntrial
    }
    return(prob.rej)
  }
  rej_null_p <- rej_null_p(pp=pp_q1[,j], cutoff_values, ntrial=1000)
  
  # Select the appropriate cutoff value based on calibration results
  selected_cutoff <- cutoff_values[which.max(rej_null_p <= 0.05)]
  selected_cutoff_list[j] <- selected_cutoff
  cat("Selected Cutoff Probability:", selected_cutoff, "\n")
}



cl<- makeCluster(23)      
registerDoParallel(cl)       

timestart <- Sys.time()

Utility <- foreach(i = 1:length(q1),.combine = 'rbind',.packages = c("LaplacesDemon","invgamma"))  %dopar% {
  #calculate parameter a and b
  para <- decide_para(x0, n0, nc, gamma, q1[i], q2, small = 0.00000000001, large = 100000000000, R = 50000)
  a <- para[1] 
  b <- para[2]
  #power
  pp1 <- ECP_normal_model(a, b, n0, x0, nc, uc=u_0, sigc=sd_0, nt, ut=u_0+0.5, sigt=sd_0, ntrial)
  rej_null <- 0
  for(l in 1:ntrial){
    if(pp1[l] > selected_cutoff_list[i]){
      rej_null <- rej_null + 1
    }
  }
  power1 <- rej_null/ntrial
  
  #error
  pp2 <- ECP_normal_model(a, b, n0, x0, nc, uc=u_0+1, sigc=sd_0, nt, ut=u_0+1, sigt=sd_0, ntrial)
  rej_null <- 0
  for(l in 1:ntrial){
    if(pp2[l] > selected_cutoff_list[i]){
      rej_null <- rej_null + 1
    }
  } 
  error1 <- rej_null/ntrial
  
  #Utility
  if (error1 > 0.05) {
    I <- 1
  } else {
    I <- 0
  }
  w2 <- 2           #Set the value of w2
  utility <- power1 - error1 - w2*(error1-0.05)*I
  return(utility)
}
stopCluster(cl)
timeend <- Sys.time()
print(timeend - timestart)
# Find the maximum value of Utility
max_value_q1 <- Utility[which.max(Utility)]


#############################################q2#######################
q1 <- 0.93  ##
q2 <- seq(0.01, 0.99, by = 0.01) 
cl<- makeCluster(23)      
registerDoParallel(cl)       
##Choose q1,q2 to maximize Utility
timestart <- Sys.time()

pp_q2 <- foreach(i = 1:length(q2),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma"))  %dopar% {
  
  #calculate parameter a and b
  para <- decide_para(x0, n0, nc, gamma, q1, q2[i], small = 0.00000000001, large = 100000000000, R = 50000)
  a <- para[1] 
  b <- para[2]
  ##probability that treatment is superior to control
  
  pp <- ECP_normal_model(a, b, n0, x0, nc, uc=u_0, sigc=sd_0, nt, ut=u_0, sigt=sd_0, ntrial)
  
  return(pp)
}
stopCluster(cl)
timeend <- Sys.time()
print(timeend - timestart)


#Calibrate Cutoff
cutoff_values <- seq(0.001, 0.999, by = 0.001)           #Set the range of cutoff values of probability that treatment is superior to control

selected_cutoff_list <- numeric(length(q2))
#function of calculating probability of rejecting null hypothesis
for (j in 1:length(q2)) {
  rej_null_p <- function(pp, cutoff_values, ntrial){
    prob.rej <- numeric(length(cutoff_values))
    for (k in 1:length(cutoff_values)) {
      cutoff <- cutoff_values[k]
      rej_null <- 0
      
      for(l in 1:ntrial){
        if(pp[l] > cutoff){
          rej_null <- rej_null + 1
        }
      }
      prob.rej[k] <- rej_null/ntrial
    }
    return(prob.rej)
  }
  rej_null_p <- rej_null_p(pp=pp_q2[,j], cutoff_values, ntrial)
  
  # Select the appropriate cutoff value based on calibration results
  selected_cutoff <- cutoff_values[which.max(rej_null_p <= 0.05)]
  selected_cutoff_list[j] <- selected_cutoff
  cat("Selected Cutoff Probability:", selected_cutoff, "\n")
}



cl<- makeCluster(23)      
registerDoParallel(cl)       

timestart <- Sys.time()

Utility <- foreach(i = 1:length(q2),.combine = 'rbind',.packages = c("LaplacesDemon","invgamma"))  %dopar% {
  #calculate parameter a and b
  para <- decide_para(x0, n0, nc, gamma, q1, q2[i], small = 0.00000000001, large = 100000000000, R = 50000)
  a <- para[1] 
  b <- para[2]
  
  #power
  pp1 <- ECP_normal_model(a, b, n0, x0, nc, uc=u_0, sigc=sd_0, nt, ut=u_0+0.5, sigt=sd_0, ntrial)
  rej_null <- 0
  for(l in 1:ntrial){
    if(pp1[l] > selected_cutoff_list[i]){
      rej_null <- rej_null + 1
    }
  }
  power1 <- rej_null/ntrial
  
  #error
  pp2 <- ECP_normal_model(a, b, n0, x0, nc, uc=u_0+1, sigc=sd_0, nt, ut=u_0+1, sigt=sd_0, ntrial)
  rej_null <- 0
  for(l in 1:ntrial){
    if(pp2[l] > selected_cutoff_list[i]){
      rej_null <- rej_null + 1
    }
  } 
  error1 <- rej_null/ntrial
  
  #Utility
  if (error1 > 0.05) {
    I <- 1
  } else {
    I <- 0
  }
  w2 <- 2           #Set the value of w2
  utility <- power1 - error1 - w2*(error1-0.05)*I
  return(utility)
}
stopCluster(cl)
timeend <- Sys.time()
print(timeend - timestart)
# Find the maximum value of Utility
max_value_q2 <- Utility[which.max(Utility)]







