###################################ECP_normal########################################
##load packages
library(dplyr)
library(LaplacesDemon)
library(invgamma)
library(foreach)
library(doParallel)
library(openxlsx)


#-----Function to decide tuning parameter a and b in elastic function g(T)=exp⁡(a+b∙log⁡(T))------------#
# Inputs:
# x0: external data
# n0: sample size for external data
# nc: sample size for control arm
# gamma: clinically highly meaningful difference
# q1: q1th percentile of congruence measure T in the homogeneous case
# q2: q2th percentile of congruence measure T in the heterogeneous case
# small: value of elastic function in the heterogeneous case
# large: value of elastic function in the homogeneous case
# R: the number of simulations
# 
# Outputs:
# a: tuning parameter in elastic function 
# b: tuning parameter in elastic function
decide_para <- function(x0, n0, nc, gamma, q1, q2, small, large, R){
  set.seed(1)
  u0 <- mean(x0)
  sig0 <- sd(x0)
  mc <- c(u0, u0 + gamma,  u0 - gamma)
  t <- matrix(NA, R, length(mc))
  for (i in 1:R) {
    for (j in 1:length(mc)) {
      xc <- rnorm(nc, mc[j], sig0)
      sp <- ((n0-1)*sig0^2 + (nc-1)*var(xc))/(n0 + nc - 2) # pooled variance
      t[i,j] <- max(n0, nc)^(-1/4)*abs(u0-mean(xc))/(sqrt(sp/n0 + sp/nc))
    }
  }
  quant1 <- quantile(t[,1], q1)
  quant2 <- quantile(t[,2], q2)
  quant3 <- quantile(t[,3], q2)
  KS_homo <- max(quant1)
  KS_hete <- min(quant2, quant3)
  b <- log(small/large)/log(KS_hete/KS_homo)
  a <- log(large)-b*log(KS_homo)
  
  return(c(a,b))
}


#-----Function to sample posterior of mean value for control arm------------#
# Inputs:
# x0: external data
# n0: sample size for external data
# xc: control data
# nc: sample size for control arm
# gt: value of elastic function at current congruence measure t between external and control data
# sim: number of simulated trial
# nburn: number of simulated trial in burn-in process
# 
# Outputs:
# muc_post: posterior of mean value for control arm
sample_poster <- function(x0, n0, xc, nc, gt, sim=20000, nburn=10000){
  muc_post <- numeric(sim-nburn)
  muc_ini <- mean(xc)
  muc <- muc_ini
  for (s in 1:sim) {
    # sample control variance from posterior
    alpha <- (1 + nc)/2   
    beta <- nc*(var(xc) + (mean(xc) - muc)^2)/2
    sig <- rinvgamma(1, alpha, beta)
    # sample mean from posterior
    D <- 1/gt+var(x0)/n0
    mu <- (nc*mean(xc)*D + sig*mean(x0))/(nc*D + sig)
    var <- sig*D/(D*nc + sig)
    muc <- rnorm(1, mu, sqrt(var)) 
    if(s > nburn){
      muc_post[s-nburn] <- muc
    }
  }
  return(list(muc_post=muc_post))
}


#-----Function to obtain the probability that treatment is superior to control-------------------#
# Inputs:
# a: calibrated parameter obtained from decide_para()
# b: calibrated parameter obtained from decide_para()
# n0: sample size for external data
# x0: external data
# nc: sample size for control arm
# uc: true mean value for control arm
# sigc: true standard deviation for control arm
# nt: sample size for treatment arm
# ut: true mean value for treatment arm
# sigt: true standard deviation for treatment arm
# ntrial: number of simulated trial
# 
# Outputs:
# pp: the list of the probability that treatment is superior to control

ECP_normal_model <- function(a, b, n0, x0, nc, uc, sigc, nt, ut, sigt, ntrial){
  pp <- numeric(ntrial)
  

  for(trial in 1:ntrial){
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rnorm(nc, uc, sigc)
    xt <- rnorm(nt, ut, sigt)
    # calculate statistic and g(t)
    sp <- ((n0-1)*var(x0) + (nc-1)*var(xc))/(n0 + nc - 2) # pooled variance
    T <- max(n0, nc)^(-1/4)*abs(mean(x0)-mean(xc))/(sqrt(sp/n0 + sp/nc))
    gt <- exp(a + b*log(T))
    if(gt < 0.000000001){gt <- 0.000000001}
    # posterior for control arm
    temp <- sample_poster(x0, n0, xc, nc, gt)
    muc <- temp$muc_post
    # posterior for treatment arm  
    t.par1 <- mean(xt)
    t.par2 <- var(xt)/nt
    mut <- rst(10000, t.par1, sqrt(t.par2), nt - 1)
    
    # probability that treatment is superior to control
    pp[trial] <- mean(mut > muc)
  }
  
  return(pp)
  
}


#function of calculating probability of rejecting null hypothesis
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
  
  rej_null_p <- cbind(cutoff_values,prob.rej)
  return(rej_null_p)
}


##########################################simulation study####################################
n0 <- 50  
u0 <- 1   
sig0 <- 1   

set.seed(8172)   ##Set random seeds to ensure reproducible results
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


#################################################################
ntrial <- 1000     #Set the number of simulations

gamma <- 1
#Determine q1 and q2 based on the grid search results
q1 <- 0.93
q2 <- 0.01
#calculate parameter a and b
para <- decide_para(x0, n0, nc, gamma, q1, q2, small = 0.00000000001, large = 100000000000, R = 50000)
a <- para[1] 
b <- para[2]

timestart <- Sys.time()
pp <- ECP_normal_model(a, b, n0, x0, nc, uc, sigc, nt, ut, sigt, ntrial)  #probability that treatment is superior to control
timeend <- Sys.time()
print(timeend - timestart)

#pp_ECP_normal_model_1000<- data.frame(pp)
#write.xlsx(pp_ECP_normal_model_1000, file = "pp_ECP_normal_model_1000.xlsx")


cutoff_values <- seq(0.001, 0.999, by = 0.001)           #Set the range of cutoff values of probability that treatment is superior to control

timestart <- Sys.time()
rej_null_p <- rej_null_p(pp, cutoff_values, ntrial)
timeend <- Sys.time()
print(timeend - timestart)
# Select the appropriate cutoff value based on calibration results
#selected_cutoff <- cutoff_values[which.min(abs(results - 0.05))]
selected_cutoff <- rej_null_p[,1][which.max(rej_null_p[,2] <= 0.05)]
cat("Selected Cutoff Probability:", selected_cutoff, "\n")


cutoff_values <- selected_cutoff          

######################Set up simulation scenarios###############
nc <- 25 
nt <- 50

sigc <-1
sigt <-1

uc <- c(1, 1  , 0.9, 1.1, 2, 2.5, 0  , -0.5)
ut <- c(1, 1.5, 1.5, 1.5, 2, 2.5, 0.5,  0   )



cl<- makeCluster(10)      
registerDoParallel(cl)     

timestart <- Sys.time()

result_pp <- foreach(i = 1:length(uc),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma","runjags"))  %dopar% {
  
  ##probability that treatment is superior to control
  
  pp <- ECP_normal_model(a, b, n0, x0, nc, uc=uc[i], sigc, nt, ut=ut[i], sigt, ntrial)
  
  return(pp)
}
stopCluster(cl)
timeend <- Sys.time()
print(timeend - timestart)

results <- numeric(length(uc))

for (j in 1:length(uc)) {
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
  results[j] <- rej_null_p(pp = result_pp[,j], cutoff_values, ntrial)
  
}

results