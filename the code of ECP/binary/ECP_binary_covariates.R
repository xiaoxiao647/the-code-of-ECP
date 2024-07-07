###################################ECP_binary_covariates########################################
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
  p <- x0/n0
  cc2 <- c(mean(cov0_2), mean(cov0_2) + gamma[3],  mean(cov0_2) - gamma[3])
  p0 <- cov0_1/n0
  if(p - gamma[1] < 0){
    mc <- c(p, p + gamma[1])
    if(p0 - gamma[2] < 0){
      cc1 <- c(p0, p0 + gamma[2])
      m <- rbind(c(mc[1], cc1[1], cc2[1]),
                 c(mc[2], cc1[2], cc2[2]),
                 c(mc[2], cc1[2], cc2[3]))
    }else if(p0 + gamma[2] > 1){
      cc1 <- c(p0, p0 - gamma[2])
      m <- rbind(c(mc[1], cc1[1], cc2[1]),
                 c(mc[2], cc1[2], cc2[2]),
                 c(mc[2], cc1[2], cc2[3]))
    }else{
      cc1 <- c(p0, p0 + gamma[2], p0 - gamma[2])
      m <- rbind(c(mc[1], cc1[1], cc2[1]),
                 c(mc[2], cc1[2], cc2[2]),
                 c(mc[2], cc1[2], cc2[3]),
                 c(mc[2], cc1[3], cc2[2]),
                 c(mc[2], cc1[3], cc2[3]))
    }
  }else if(p + gamma[1] > 1){
    mc <- c(p, p - gamma[1])
    if(p0 - gamma[2] < 0){
      cc1 <- c(p0, p0 + gamma[2])
      m <- rbind(c(mc[1], cc1[1], cc2[1]),
                 c(mc[2], cc1[2], cc2[2]),
                 c(mc[2], cc1[2], cc2[3]))
    }else if(p0 + gamma[2] > 1){
      cc1 <- c(p0, p0 - gamma[2])
      m <- rbind(c(mc[1], cc1[1], cc2[1]),
                 c(mc[2], cc1[2], cc2[2]),
                 c(mc[2], cc1[2], cc2[3]))
    }else{
      cc1 <- c(p0, p0 + gamma[2], p0 - gamma[2])
      m <- rbind(c(mc[1], cc1[1], cc2[1]),
                 c(mc[2], cc1[2], cc2[2]),
                 c(mc[2], cc1[2], cc2[3]),
                 c(mc[2], cc1[3], cc2[2]),
                 c(mc[2], cc1[3], cc2[3]))
    }
  }else{
    mc <- c(p, p + gamma[1], p - gamma[1])
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
  }
  t <- matrix(NA, R, dim(m)[1])
  for (i in 1:R) {
    for (j in 1:dim(m)[1]) {
      out <- rbinom(1, nc, m[j, 1])
      cov1 <- rbinom(1, nc, m[j, 2])
      cov2 <- rnorm(nc, m[j, 3], sd(cov0_2))
      phat <- (x0 + out)/(n0 + nc)
      O <- cbind(c(x0, out), c(n0-x0, nc-out))
      E <- cbind(c(n0, nc) * phat, c(n0, nc) * (1 - phat))
      chi <- sum((O - E)^2/E)
      p1 <- 1-pchisq(max(n0, nc)^(-1/4)*chi, 1)
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

#Defeine JAGS model
ECP_binary_modle =  "
model {
  #historical likelihood
  x0 ~ dbin(p0,n0) 
  
  #current likelihood
  xc ~ dbin(pc,nc)
  
  logit_pc ~ dnorm(log(p0/(1-p0)),gt)
  pc <- exp(logit_pc)/(1+exp(logit_pc))
  
  p0 ~ dbeta(1/2,1/2)

} 
"
writeLines(ECP_binary_modle,con = "ECP_binary_modle.txt")


ECP_cov_binary <- function(n0, x0, cov0_1, cov0_2, nc, pc, pc_1, cu, csig, nt, pt, a, b, t_alpha, t_beta, ntrial){
  pp <- numeric(ntrial)
  
  for(trial in 1:ntrial){
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rbinom(1, nc, pc)
    xt <- rbinom(1, nt, pt)
    cov1 <- rbinom(1, nc, pc_1)
    cov2 <- rnorm(nc, cu, csig)
    # calculate statistic and g(t)
    # obtain p value for outcome
    phat <- (x0 + xc)/(n0 + nc)
    O <- cbind(c(x0, xc), c(n0-x0, nc-xc))
    E <- cbind(c(n0, nc) * phat, c(n0, nc) * (1 - phat))
    chi <- sum((O - E)^2/E)
    p1 <- 1-pchisq(max(n0, nc)^(-1/4)*chi, 1)
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
    if(gt > 10000){gt <- 10000}
    
    #Data packing
    dat <- dump.format(list(xc=xc, nc=nc, n0=n0, x0=x0, gt=gt))
    # posterior for control arm
    ##Initialize chains
    inits1 = dump.format(list( .RNG.name="base::Super-Duper", .RNG.seed=99999))
    inits2 = dump.format(list( .RNG.name="base::Wichmann-Hill", .RNG.seed=12345))
    inits3 = dump.format(list( .RNG.name="base::Mersenne-Twister", .RNG.seed=1802))
    
    ##Tell JAGS which latent variables to monitor
    monitor = c("pc","p0")
    
    ##Run JAGS
    results <- run.jags(model = "ECP_binary_modle.txt",
                        monitor = monitor,
                        data = dat,
                        n.chains = 3,
                        inits = c(inits1,inits2,inits3),
                        plots = FALSE,
                        burnin = 10000,
                        sample = 10000,
                        thin = 1)
    chains = rbind(results$mcmc[[1]],results$mcmc[[2]],results$mcmc[[3]])
    pc_post = sample(chains[,"pc"],10000)
    # posterior for treatment arm  
    pt_post <- rbeta(10000, t_alpha + xt, t_beta + nt - xt)
    
    # probability that treatment is superior to control
    pp[trial] <- mean(pt_post > pc_post)
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
n0 <- 100  
p0 <- 0.4   

p_cov0_1 <- 0.4  #True probability value of covariate 1 in external data

u_cov0_2 <- 4    #True mean of covariate 2 in external data
sig_cov0_2 <-1   #True variance of covariate 2 in external data

set.seed(121)
x0 <- rbinom(1, n0, p0)

set.seed(121)
cov0_1 <- rbinom(1, n0, p_cov0_1)

set.seed(20)
cov0_2 <- rnorm(n0, u_cov0_2, sig_cov0_2)


## control and treatment
nc <- 40
nt <- 80

pc <- x0/n0
pt <- x0/n0

pc_1 <- cov0_1/n0
cu   <- mean(cov0_2)
csig <- sd(cov0_2)

#################################################################
ntrial <- 1000     #Set the number of simulations
t_alpha <- 1/2
t_beta  <- 1/2
#calculate parameter a and b
para <- decide_para(n0, x0, cov0_1, cov0_2, nc, gamma=c(0.2, 0.2, 1.5), q1=0.999, q2=0.001, small = 0.00000000001, large = 100000000000, R = 50000)
a <- para[1] 
b <- para[2]

timestart <- Sys.time()
pp <- ECP_cov_binary(n0, x0, cov0_1, cov0_2, nc, pc, pc_1, cu, csig, nt, pt, a, b, t_alpha, t_beta, ntrial)  #probability that treatment is superior to control
timeend <- Sys.time()
print(timeend - timestart)

pp_ECP_cov_binary_0.999_0.001_1000<- data.frame(pp)
write.xlsx(pp_ECP_cov_binary_0.999_0.001_1000, file = "pp_ECP_cov_binary_0.999_0.001_1000.xlsx")


cutoff_values <- seq(0.001, 0.999, by = 0.001)            #Set the range of cutoff values of probability that treatment is superior to control

timestart <- Sys.time()
rej_null_p <- rej_null_p(pp, cutoff_values, ntrial)
timeend <- Sys.time()
print(timeend - timestart)
# Select the appropriate cutoff value based on calibration results
#selected_cutoff <- cutoff_values[which.min(abs(results - 0.05))]
selected_cutoff <- rej_null_p[,1][which.max(rej_null_p[,2] <= 0.05)]
cat("Selected Cutoff Probability:", selected_cutoff, "\n")

cutoff_values <- 0.92            

######################Set up simulation scenarios###############
nc <- 40 
nt <- 80

pc <- c(0.4, 0.4, 0.35, 0.45, 0.6, 0.65, 0.2, 0.2)
pt <- c(0.4, 0.6, 0.55, 0.6, 0.6, 0.65, 0.35, 0.4)

pc_1 <- c(0.4,0.4,0.38,0.42,0.6,0.66,0.2,0.15)

cu <- c(4,4,3.9,4.1,5.5,6,2.5,2)
csig <- 1

ntrial <- 1000  

cl<- makeCluster(10)      
registerDoParallel(cl)       

timestart <- Sys.time()

result_pp <- foreach(i = 1:length(pc),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma","runjags","metap"))  %dopar% {
  
  ##probability that treatment is superior to control
  
  pp <- ECP_cov_binary(n0, x0, cov0_1, cov0_2, nc, pc=pc[i], pc_1=pc_1[i], cu=cu[i], csig, nt, pt=pt[i], a, b, t_alpha, t_beta, ntrial)
  
  return(pp)
}
stopCluster(cl)
timeend <- Sys.time()
print(timeend - timestart)

results <- numeric(length(pc))

for (j in 1:length(pc)) {
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

