###################################ECP_normal_covariates########################################
##load packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")


library(dplyr)
library(LaplacesDemon)
library(invgamma)
library(metap)
library(foreach)
library(doParallel)
library(openxlsx)


#-----Function to decide tuning parameter a and b in elastic function g(T)=exp⁡(a+b∙log⁡(T))------------#
# Inputs:
# n0: sample size for external data
# x0: outcome data for external arm
# cov0_1: covariate 1 for external arm
# cov0_2: covariate 2 for external arm
# nc: sample size for control arm
# gamma: clinically highly meaningful differences for outcome and covariates
# q1: q1th percentile of congruence measure T in the homogeneous case
# q2: q2th percentile of congruence measure T in the heterogeneous case
# small: value of elastic function in the heterogeneous case
# large: value of elastic function in the homogeneous case
# R: the number of simulations
# 
# Outputs:
# a: tuning parameter in elastic function 
# b: tuning parameter in elastic function
decide_para <- function(n0, x0, cov0_1, cov0_2, nc, gamma, q1, q2, small, large, R=50000){
  set.seed(1)
  mc <- c(mean(x0), mean(x0) + gamma[1],  mean(x0) - gamma[1])
  cc2 <- c(mean(cov0_2), mean(cov0_2) + gamma[3],  mean(cov0_2) - gamma[3])
  p0 <- cov0_1/n0
  if(p0 - gamma[2] < 0){
    cc1 <- c(p0, p0 + gamma[2])
    m <- rbind(c(mc[1], cc1[1], cc2[1]),
               c(mc[2], cc1[2], cc2[2]),
               c(mc[2], cc1[2], cc2[3]),
               c(mc[3], cc1[2], cc2[2]),
               c(mc[3], cc1[2], cc2[3]))
  }else if(p0 + gamma[2] > 1){
    cc1 <- c(p0, p0 - gamma[2])
    m <- rbind(c(mc[1], cc1[1], cc2[1]),
               c(mc[2], cc1[2], cc2[2]),
               c(mc[2], cc1[2], cc2[3]),
               c(mc[3], cc1[2], cc2[2]),
               c(mc[3], cc1[2], cc2[3]))
  }else{
    cc1 <- c(p0, p0 + gamma[2], p0 - gamma[2])
    m <- rbind(c(mc[1], cc1[1], cc2[1]),
               c(mc[2], cc1[2], cc2[2]),
               c(mc[2], cc1[2], cc2[3]),
               c(mc[2], cc1[3], cc2[2]),
               c(mc[3], cc1[2], cc2[2]),
               c(mc[2], cc1[3], cc2[3]),
               c(mc[3], cc1[2], cc2[3]),
               c(mc[3], cc1[3], cc2[2]),
               c(mc[3], cc1[3], cc2[3]))
  }
  t <- matrix(NA, R, dim(m)[1])
  for (i in 1:R) {
    for (j in 1:dim(m)[1]) {
      out <- rnorm(nc, m[j, 1], sd(x0))
      cov1 <- rbinom(1, nc, m[j, 2])
      cov2 <- rnorm(nc, m[j, 3], sd(cov0_2))
      s <- sqrt(((n0-1)*var(x0) + (nc-1)*var(out))/(n0+nc-2))
      tstat <- abs((mean(x0)-mean(out))/(s*sqrt(1/n0+1/nc)))
      p1 <-2*(1-pt(max(n0, nc)^(-1/4)*tstat, nc+n0-2))
      phat <- (cov0_1 + cov1)/(n0 + nc)
      O <- cbind(c(cov0_1, cov1), c(n0-cov0_1, nc-cov1))
      E <- cbind(c(n0, nc) * phat, c(n0, nc) * (1 - phat))
      chi <- sum((O - E)^2/E)
      p2 <- 1-pchisq(max(n0, nc)^(-1/4)*chi, 1)
      s <- sqrt(((n0-1)*var(cov0_2) + (nc-1)*var(cov2))/(n0+nc-2))
      tstat <- abs((mean(cov0_2)-mean(cov2))/(s*sqrt(1/n0+1/nc)))
      p3 <-2*(1-pt(max(n0, nc)^(-1/4)*tstat, nc+n0-2))
      pval <- c(p1, p2, p3)
      t[i,j] <- -log(sumlog(pval)$p)
    }
  }
  quant1 <- quantile(t[,1], q1)
  quant2 <- sapply(2:dim(t)[2], function(x) quantile(t[ ,x], q2))
  KS_homo <- max(quant1)
  KS_hete <- min(quant2)
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
    D <- var(x0)/(n0*gt)
    mu <- (nc*mean(xc)*D + sig*mean(x0))/(nc*D + sig)
    var <- sig*D/(D*nc + sig)
    muc <- rnorm(1, mu, sqrt(var)) 
    if(s > nburn){
      muc_post[s-nburn] <- muc
    }
  }
  return(muc_post=muc_post)
}


#-----Function to obtain the probability that treatment is superior to control-------------------#
# Inputs:
# n0: sample size for external data
# x0: outcome data for external arm
# cov0_1: covariate 1 for external arm
# cov0_2: covariate 2 for external arm
# nc: sample size for control arm
# uc: true mean value of outcome for control arm
# sigc: true standard deviation of outcome for control arm
# pc_1: true rate of covariate 1 for control arm
# cu: true mean value of covariate 2 for control arm
# csig: true standard deviarion of covariate 2 for control arm
# nt: sample size for treatment arm
# ut: true mean value of outcome for treatment arm
# sigt: true standard deviation of outcome for treatment arm
# a: calibrated parameter obtained from decide_para()
# b: calibrated parameter obtained from decide_para()
# ntrial: number of simulated trial
# 
# Outputs:
# pp: the list of the probability that treatment is superior to control

ECP_cov_normal <- function(n0, x0, cov0_1, cov0_2, nc, uc, sigc, pc_1, cu, csig, nt, ut, sigt, a, b, ntrial){
  pp <- numeric(ntrial)
  
  for(trial in 1:ntrial){
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rnorm(nc, uc, sigc)
    xt <- rnorm(nt, ut, sigt)
    cov1 <- rbinom(1, nc, pc_1)
    cov2 <- rnorm(nc, cu, csig)
    # calculate statistic and g(t)
    # obtain p value for outcome
    s <- sqrt(((n0-1)*var(x0) + (nc-1)*var(xc))/(n0+nc-2))
    tstat <- abs((mean(x0)-mean(xc))/(s*sqrt(1/n0+1/nc)))
    p1 <-2*(1-pt(max(n0, nc)^(-1/4)*tstat, nc+n0-2))
    # obtain p value for covariate 1
    phat <- (cov0_1 + cov1)/(n0 + nc)
    O <- cbind(c(cov0_1, cov1), c(n0-cov0_1, nc-cov1))
    E <- cbind(c(n0, nc) * phat, c(n0, nc) * (1 - phat))
    chi <- sum((O - E)^2/E)
    p2 <- 1-pchisq(max(n0, nc)^(-1/4)*chi, 1)
    # obtain p value for covariate 2
    s <- sqrt(((n0-1)*var(cov0_2) + (nc-1)*var(cov2))/(n0+nc-2))
    tstat <- abs((mean(cov0_2)-mean(cov2))/(s*sqrt(1/n0+1/nc)))
    p3 <-2*(1-pt(max(n0, nc)^(-1/4)*tstat, nc+n0-2))
    pval <- c(p1, p2, p3)
    T <- -log(sumlog(pval)$p)
    gt <- exp(a + b*log(T))
    if(gt < 0.00001){gt <- 0.00001}
    # posterior for control arm
    muc <- sample_poster(x0, n0, xc, nc, gt)
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

p_cov0_1 <- 0.4  #True probability value of covariate 1 in external data

u_cov0_2 <- 4    #True mean of covariate 2 in external data
sig_cov0_2 <-1   #True variance of covariate 2 in external data

set.seed(8172)
x0 <- rnorm(n0, u0, sig0)

set.seed(113)
cov0_1 <- rbinom(1, n0, p_cov0_1)

set.seed(4310)
cov0_2 <- rnorm(n0, u_cov0_2, sig_cov0_2)


## control and treatment
nc <- 25
nt <- 50

uc <- mean(x0)
ut <- mean(x0)
sigc <- sd(x0)
sigt <- sd(x0)

pc_1 <- cov0_1/n0
cu   <- mean(cov0_2)
csig <- sd(cov0_2)


#################################################################
ntrial <- 1000     #Set the number of simulations

#calculate parameter a and b
para <- decide_para(n0, x0, cov0_1, cov0_2, nc, gamma=c(1, 0.2, 1.5), q1=0.999, q2=0.001, small = 0.00000000001, large = 100000000000, R = 50000)
a <- para[1] 
b <- para[2]

timestart <- Sys.time()
pp <- ECP_cov_normal(n0, x0, cov0_1, cov0_2, nc, uc, sigc, pc_1, cu, csig, nt, ut, sigt, a, b, ntrial)  #probability that treatment is superior to control
timeend <- Sys.time()
print(timeend - timestart)

#pp_ECP_cov_normal_0.999_0.001_1000<- data.frame(pp)
#write.xlsx(pp_ECP_cov_normal_0.999_0.001_1000, file = "pp_ECP_cov_normal_0.999_0.001_1000.xlsx")


cutoff_values <- seq(0.001, 0.999, by = 0.001)            #Set the range of cutoff values of probability that treatment is superior to control


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
sigc <- 1
sigt <- 1

uc <- c(1, 1, 0.9, 1.1, 2, 2.5, 0, -0.5)
ut <- c(1, 1.5, 1.5, 1.5, 2, 2.5, 0.5, 0)

pc_1 <- c(0.4,0.4,0.38,0.42,0.6,0.66,0.2,0.15)

cu <- c(4,4,3.9,4.1,5.5,6,2.5,2)
csig <- 1

ntrial <- 1000  

cl<- makeCluster(10)      
registerDoParallel(cl)   

timestart <- Sys.time()

result_pp <- foreach(i = 1:length(muc),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma","runjags","metap"))  %dopar% {

  ##probability that treatment is superior to control
  
  pp <- ECP_cov_normal(n0, x0, cov0_1, cov0_2, nc, uc=uc[i], sigc, pc_1=pc_1[i], cu=cu[i], csig, nt, ut=ut[i], sigt, a, b, ntrial)
  
  return(pp)
}
stopCluster(cl)
timeend <- Sys.time()
print(timeend - timestart)

results <- numeric(length(uc))

for (j in 1:length(muc)) {
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

