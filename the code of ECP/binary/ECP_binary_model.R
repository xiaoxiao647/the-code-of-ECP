###################################ECP_binary########################################
##load packages
library(LaplacesDemon)
library(runjags)
library(openxlsx)

library(dplyr)
library(invgamma)
library("foreach")
library("doParallel")

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
  set.seed(2)
  p0 <- x0/n0
  if(p0-gamma<0){
    p <- c(p0, p0 + gamma)
  }else if(p0 + gamma>1){
    p <- c(p0, p0 - gamma)
  }else{
    p <- c(p0, p0 - gamma, p0 + gamma)
  }
  K <- matrix(NA, R, length(p))
  for (i in 1:R) {
    for (j in 1:length(p)) {
      y <- rbinom(1, nc, p[j])
      phat <- (x0 + y)/(n0 + nc)
      obs <- cbind(c(x0, y), c(n0-x0, nc-y))
      exc <- cbind(c(n0, nc)*phat, c(n0, nc)*(1-phat))
      Ka <- max(n0, nc)^(-1/4)*sum((obs-exc)^2/exc)
      K[i,j] <- Ka
    }
  }
  if(length(p)==2){
    quant1 <- quantile(K[,1], probs = q1)
    quant2 <- quantile(K[,2], probs = q2)
    K_homo <- quant1
    K_hete <- quant2
  }else{
    quant1 <- quantile(K[,1], probs = q1)
    quant2 <- quantile(K[,2], probs = q2)
    quant3 <- quantile(K[,3], probs = q2)
    K_homo <- max(quant1)
    K_hete <- min(quant2, quant3)
  }
  b <- log(small/large)/log(K_hete/K_homo)
  a <- log(large)-b*log(K_homo)
  
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



################################ECP_binary_modle#################################
ECP_binary_modle <- function(a, b, n0, x0, nc, pc, nt, pt, t_alpha, t_beta, ntrial){
  pp <- numeric(ntrial)
  
  for(trial in 1:ntrial){
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rbinom(1, nc, pc)
    xt <- rbinom(1, nt, pt)
    
    #calculate g(t)
    phat <- (x0 + xc)/(n0 + nc)
    O <- cbind(c(x0, xc), c(n0-x0, nc-xc))
    E <- cbind(c(n0, nc) * phat, c(n0, nc) * (1 - phat))
    T <- max(n0, nc)^(-1/4)*sum((O - E)^2/E)
    gt <- exp(a + b*log(T))
    if(gt < 0.00001) gt <- 0.00001
    if(gt > 10000) gt <- 10000
    
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

set.seed(121)     ##Set random seeds to ensure reproducible results
x0 <- rbinom(1, n0, p0)

## control and treatment
nc <- 40
nt <- 80

pc <- x0/n0
pt <- x0/n0

#################################################################
ntrial <- 1000    #Set the number of simulations
t_alpha <- 1/2
t_beta  <- 1/2

gamma <- 0.2      
#Determine q1 and q2 based on the grid search results
q1 <- 0.91          
q2 <- 0.22          
## obtain calibrated a and b in elastic function g(T)=exp⁡(a+b∙log⁡(T))
para <- decide_para(x0, n0, nc, gamma, q1, q2, small = 0.00000000001, large = 100000000000, R = 50000)
a <- para[1] 
b <- para[2]


timestart <- Sys.time()
pp <- ECP_binary_modle(a, b, n0, x0, nc, pc, nt, pt, t_alpha, t_beta, ntrial)  #probability that treatment is superior to control
timeend <- Sys.time()
print(timeend - timestart)

pp_ECP_binary_1000<- data.frame(pp)
write.xlsx(pp_ECP_binary_1000, file = "pp_ECP_binary_1000.xlsx")


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
nc <- 40 
nt <- 80

pc <- c(0.4, 0.4, 0.35, 0.45, 0.6, 0.65, 0.2, 0.2)
pt <- c(0.4, 0.6, 0.55, 0.6 , 0.6, 0.65, 0.35, 0.4)


## obtain calibrated a and b in elastic function g(T)=exp⁡(a+b∙log⁡(T))
para <- decide_para(x0, n0, nc, gamma, q1, q2, small = 0.00000000001, large = 100000000000, R = 50000)
a <- para[1] 
b <- para[2]

cl<- makeCluster(10)      
registerDoParallel(cl)      

timestart <- Sys.time()

result_pp <- foreach(i = 1:length(pc),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma","runjags"))  %dopar% {
  
  ##probability that treatment is superior to control
  
  pp <- ECP_binary_modle(a, b, n0, x0, nc, pc=pc[i], nt, pt=pt[i], t_alpha, t_beta, ntrial)
  
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