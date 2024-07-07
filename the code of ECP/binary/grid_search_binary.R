##########################################grid_search_binary#################################################

#Generate simulated external data
n0 <-100  # Sample size of the external control group
p0 <- 0.4   # True probability value for the external control group

## Set random seeds to ensure reproducible results
set.seed(121)
x0 <- rbinom(1, n0, p0)

#Control and treatment

nc <- 40
nt <- 80

pc <- x0/n0
pt <- x0/n0


ntrial <- 100     #Set the number of simulations
t_alpha <- 1/2
t_beta  <- 1/2

gamma <- 0.2       #Set the value of gamma
#############################################q1###########################
q1 <- seq(0.01, 0.99, by = 0.01) ##
q2 <-  0.1

cl<- makeCluster(23)      
registerDoParallel(cl)   

##Choose q1,q2 to maximize Utility
timestart <- Sys.time()

pp_q1 <- foreach(i = 1:length(q1),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma","runjags"))  %dopar% {
  
  #calculate parameter a and b
  para <- decide_para(x0, n0, nc, gamma, q1[i], q2, small = 0.00000000001, large = 100000000000, R = 50000)
  a <- para[1] 
  b <- para[2]
  ##probability that treatment is superior to control
  
  pp <- ECP_binary_modle(a, b, n0, x0, nc, pc, nt, pt, t_alpha, t_beta, ntrial)  #probability that treatment is superior to control
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
  rej_null_p <- rej_null_p(pp=pp_q1[,j], cutoff_values, ntrial)
  
  # Select the appropriate cutoff value based on calibration results
  selected_cutoff <- cutoff_values[which.max(rej_null_p <= 0.05)]
  selected_cutoff_list[j] <- selected_cutoff
  cat("Selected Cutoff Probability:", selected_cutoff, "\n")
}



cl<- makeCluster(23)      
registerDoParallel(cl)       


timestart <- Sys.time()

Utility <- foreach(i = 1:length(q1),.combine = 'rbind',.packages = c("LaplacesDemon","invgamma","runjags"))  %dopar% {
  #calculate parameter a and b
  para <- decide_para(x0, n0, nc, gamma, q1[i], q2, small = 0.00000000001, large = 100000000000, R = 50000)
  a <- para[1] 
  b <- para[2]
  #power
  pp1 <-  ECP_binary_modle(a, b, n0, x0, nc, pc, nt, pt+gamma, t_alpha, t_beta, ntrial)  #probability that treatment is superior to control
  rej_null <- 0
  for(l in 1:ntrial){
    if(pp1[l] > selected_cutoff_list[i]){
      rej_null <- rej_null + 1
    }
  }
  power1 <- rej_null/ntrial
  
  #error
  pp2 <-  ECP_binary_modle(a, b, n0, x0, nc, pc+gamma, nt, pt+gamma, t_alpha, t_beta, ntrial)   #probability that treatment is superior to control
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
q1 <- 0.99  ##
q2 <- seq(0.01, 0.99, by = 0.01) 
cl<- makeCluster(23)      
registerDoParallel(cl)       
##Choose q1,q2 to maximize Utility
timestart <- Sys.time()

pp_q2 <- foreach(i = 1:length(q2),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma","runjags"))  %dopar% {
  
  #calculate parameter a and b
  para <- decide_para(x0, n0, nc, gamma, q1, q2[i], small = 0.00000000001, large = 100000000000, R = 50000)
  a <- para[1] 
  b <- para[2]
  ##probability that treatment is superior to control
  
  pp <- ECP_binary_modle(a, b, n0, x0, nc, pc, nt, pt, t_alpha, t_beta, ntrial)  #probability that treatment is superior to control
  
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
  rej_null_p <- rej_null_p(pp=pp_q2[,j], cutoff_values, ntrial=100)
  
  # Select the appropriate cutoff value based on calibration results
  selected_cutoff <- cutoff_values[which.max(rej_null_p <= 0.05)]
  selected_cutoff_list[j] <- selected_cutoff
  cat("Selected Cutoff Probability:", selected_cutoff, "\n")
}



cl<- makeCluster(23)      
registerDoParallel(cl)      

timestart <- Sys.time()

Utility <- foreach(i = 1:length(q2),.combine = 'rbind',.packages = c("LaplacesDemon","invgamma","runjags"))  %dopar% {
  #calculate parameter a and b
  para <- decide_para(x0, n0, nc, gamma, q1, q2[i], small = 0.00000000001, large = 100000000000, R = 50000)
  a <- para[1] 
  b <- para[2]
  
  #power
  pp1 <- ECP_binary_modle(a, b, n0, x0, nc, pc, nt, pt+gamma, t_alpha, t_beta, ntrial)  #probability that treatment is superior to control
  rej_null <- 0
  for(l in 1:ntrial){
    if(pp1[l] > selected_cutoff_list[i]){
      rej_null <- rej_null + 1
    }
  }
  power1 <- rej_null/ntrial
  
  #error
  pp2 <- ECP_binary_modle(a, b, n0, x0, nc, pc+gamma, nt, pt+gamma, t_alpha, t_beta, ntrial)  #probability that treatment is superior to control
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

