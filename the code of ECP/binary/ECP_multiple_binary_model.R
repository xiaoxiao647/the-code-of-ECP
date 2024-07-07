###################################ECP_multiple_binary########################################
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
ECP_multiple_p0_binary_modle =  "
model {
  
  #historical likelihood
  x0  ~ dbin(p0,n0)
  
  p0  ~ dbeta(1/2, 1/2)
  mu0 <- log(p0/(1-p0))

} 
"
writeLines(ECP_multiple_p0_binary_modle,con = "ECP_multiple_p0_binary_modle.txt")





#Defeine JAGS model
ECP_multiple_binary_modle =  "
model {
  
  #current likelihood
  xc  ~ dbin(pc,nc)
  
  logit_pc ~ dnorm(mu_sum, 1/sig_sum)
  pc <- exp(logit_pc)/(1+exp(logit_pc))

} 
"
writeLines(ECP_multiple_binary_modle,con = "ECP_multiple_binary_modle.txt")



################################ECP_binary_modle#################################
ECP_multiple_binary_modle <- function(a01, b01, a02, b02, a03, b03, a04, b04, n01, x01, n02, x02, n03, x03, n04, x04, mu01, sig01, mu02, sig02, mu03, sig03, mu04, sig04, nc, pc, nt, pt, t_alpha, t_beta, ntrial){
  pp <- numeric(ntrial)
  
  for(trial in 1:ntrial){
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rbinom(1, nc, pc)
    xt <- rbinom(1, nt, pt)
    
    #calculate g(t)
    phat01 <- (x01 + xc)/(n01 + nc)
    O01 <- cbind(c(x01, xc), c(n01-x01, nc-xc))
    E01 <- cbind(c(n01, nc) * phat01, c(n01, nc) * (1 - phat01))
    T01 <- max(n01, nc)^(-1/4)*sum((O01 - E01)^2/E01)
    gt01 <- exp(a01 + b01*log(T01))
    if(gt01 < 0.000000001) gt01 <- 0.000000001
    if(gt01 > 10000) gt01 <- 10000
    
    phat02 <- (x02 + xc)/(n02 + nc)
    O02 <- cbind(c(x02, xc), c(n02-x02, nc-xc))
    E02 <- cbind(c(n02, nc) * phat02, c(n02, nc) * (1 - phat02))
    T02 <- max(n02, nc)^(-1/4)*sum((O02 - E02)^2/E02)
    gt02 <- exp(a02 + b02*log(T02))
    if(gt02 < 0.000000001) gt02 <- 0.000000001
    if(gt02 > 10000) gt02 <- 10000
    
    phat03 <- (x03 + xc)/(n03 + nc)
    O03 <- cbind(c(x03, xc), c(n03-x03, nc-xc))
    E03 <- cbind(c(n03, nc) * phat03, c(n03, nc) * (1 - phat03))
    T03 <- max(n03, nc)^(-1/4)*sum((O03 - E03)^2/E03)
    gt03 <- exp(a03 + b03*log(T03))
    if(gt03 < 0.000000001) gt03 <- 0.000000001
    if(gt03 > 10000) gt03 <- 10000
    
    phat04 <- (x04 + xc)/(n04 + nc)
    O04 <- cbind(c(x04, xc), c(n04-x04, nc-xc))
    E04 <- cbind(c(n04, nc) * phat04, c(n04, nc) * (1 - phat04))
    T04 <- max(n04, nc)^(-1/4)*sum((O04 - E04)^2/E04)
    gt04 <- exp(a04 + b04*log(T04))
    if(gt04 < 0.000000001) gt04 <- 0.000000001
    if(gt04 > 10000) gt04 <- 10000
    
    tau01 <- 1/(1/gt01+sig01/n01)
    tau02 <- 1/(1/gt02+sig02/n02)
    tau03 <- 1/(1/gt03+sig03/n03)
    tau04 <- 1/(1/gt04+sig04/n04)
    sig_sum <- 1/(tau01 + tau02 + tau03 + tau04)
    mu_sum  <- (tau01*mu01 + tau02*mu02 + tau03*mu03 + tau04*mu04) * sig_sum
    
    
    # posterior for control arm
    #Data packing
    dat <- dump.format(list(xc=xc, nc=nc, sig_sum=sig_sum, mu_sum=mu_sum))
    ##Initialize chains
    inits1 = dump.format(list( .RNG.name="base::Super-Duper", .RNG.seed=99999))
    inits2 = dump.format(list( .RNG.name="base::Wichmann-Hill", .RNG.seed=12345))
    inits3 = dump.format(list( .RNG.name="base::Mersenne-Twister", .RNG.seed=1802))
    
    ##Tell JAGS which latent variables to monitor
    monitor = c("pc")
    
    ##Run JAGS
    results <- run.jags(model = "ECP_multiple_binary_modle.txt",
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

## genertae four external data
set.seed(1)
p0 <- 0.4
u0 <- log(p0/(1-p0))
mu <- rnorm(4, u0, 0.1)
p <- exp(mu)/(1+exp(mu))

set.seed(14767)
# generate external data1
n01 <- 80
x01 <- rbinom(1, n01, p[1])

set.seed(13501)
# generate external data2
n02 <- 100
x02 <- rbinom(1, n02, p[2])

set.seed(4640)
# generate external data3
n03 <- 90
x03 <- rbinom(1, n03, p[3])

set.seed(296)
# generate external data4
n04 <- 110
x04 <- rbinom(1, n04, p[4])


#Data packing
dat <- dump.format(list(x0=x01, n0=n01))
##Initialize chains
inits1 = dump.format(list( .RNG.name="base::Super-Duper", .RNG.seed=99999))
inits2 = dump.format(list( .RNG.name="base::Wichmann-Hill", .RNG.seed=12345))
inits3 = dump.format(list( .RNG.name="base::Mersenne-Twister", .RNG.seed=1802))

##Tell JAGS which latent variables to monitor
monitor = c("mu0")

##Run JAGS
results01 <- run.jags(model = "ECP_multiple_p0_binary_modle.txt",
                    monitor = monitor,
                    data = dat,
                    n.chains = 3,
                    inits = c(inits1,inits2,inits3),
                    plots = FALSE,
                    burnin = 10000,
                    sample = 10000,
                    thin = 1)
mu01 <- results01$summary[1]$statistics[1,2]
sig01 <- results01$summary[1]$statistics[1,2]^2


#Data packing
dat <- dump.format(list(x0=x02, n0=n02))
##Initialize chains
inits1 = dump.format(list( .RNG.name="base::Super-Duper", .RNG.seed=99999))
inits2 = dump.format(list( .RNG.name="base::Wichmann-Hill", .RNG.seed=12345))
inits3 = dump.format(list( .RNG.name="base::Mersenne-Twister", .RNG.seed=1802))

##Tell JAGS which latent variables to monitor
monitor = c("mu0")

##Run JAGS
results02 <- run.jags(model = "ECP_multiple_p0_binary_modle.txt",
                      monitor = monitor,
                      data = dat,
                      n.chains = 3,
                      inits = c(inits1,inits2,inits3),
                      plots = FALSE,
                      burnin = 10000,
                      sample = 10000,
                      thin = 1)
mu02 <- results02$summary[1]$statistics[1,2]
sig02 <- results02$summary[1]$statistics[1,2]^2

#Data packing
dat <- dump.format(list(x0=x03, n0=n03))
##Initialize chains
inits1 = dump.format(list( .RNG.name="base::Super-Duper", .RNG.seed=99999))
inits2 = dump.format(list( .RNG.name="base::Wichmann-Hill", .RNG.seed=12345))
inits3 = dump.format(list( .RNG.name="base::Mersenne-Twister", .RNG.seed=1802))

##Tell JAGS which latent variables to monitor
monitor = c("mu0")

##Run JAGS
results03 <- run.jags(model = "ECP_multiple_p0_binary_modle.txt",
                      monitor = monitor,
                      data = dat,
                      n.chains = 3,
                      inits = c(inits1,inits2,inits3),
                      plots = FALSE,
                      burnin = 10000,
                      sample = 10000,
                      thin = 1)
mu03 <- results03$summary[1]$statistics[1,2]
sig03 <- results03$summary[1]$statistics[1,2]^2

#Data packing
dat <- dump.format(list(x0=x04, n0=n04))
##Initialize chains
inits1 = dump.format(list( .RNG.name="base::Super-Duper", .RNG.seed=99999))
inits2 = dump.format(list( .RNG.name="base::Wichmann-Hill", .RNG.seed=12345))
inits3 = dump.format(list( .RNG.name="base::Mersenne-Twister", .RNG.seed=1802))

##Tell JAGS which latent variables to monitor
monitor = c("mu0")

##Run JAGS
results04 <- run.jags(model = "ECP_multiple_p0_binary_modle.txt",
                      monitor = monitor,
                      data = dat,
                      n.chains = 3,
                      inits = c(inits1,inits2,inits3),
                      plots = FALSE,
                      burnin = 10000,
                      sample = 10000,
                      thin = 1)
mu04 <- results04$summary[1]$statistics[1,2]
sig04 <- results04$summary[1]$statistics[1,2]^2

## control and treatment
nc <- 40
nt <- 80

pc <- (x01/n01+x02/n02+x03/n03+x04/n04)/4
pt <- (x01/n01+x02/n02+x03/n03+x04/n04)/4



gamma <- 0.2       
#Determine q1 and q2 based on the grid search results
q1 <- 0.88          
q2 <- 0.22          
## obtain calibrated a and b in elastic function g(T)=exp⁡(a+b∙log⁡(T))
para01 <- decide_para(x0 = x01, n0 = n01, nc, gamma, q1, q2, small = 0.00000000001, large = 100000000000, R = 50000)
para02 <- decide_para(x0 = x02, n0 = n02, nc, gamma, q1, q2, small = 0.00000000001, large = 100000000000, R = 50000)
para03 <- decide_para(x0 = x03, n0 = n03, nc, gamma, q1, q2, small = 0.00000000001, large = 100000000000, R = 50000)
para04 <- decide_para(x0 = x04, n0 = n04, nc, gamma, q1, q2, small = 0.00000000001, large = 100000000000, R = 50000)

a01 <- para01[1] 
b01 <- para01[2]
a02 <- para02[1] 
b02 <- para02[2]
a03 <- para03[1] 
b03 <- para03[2]
a04 <- para04[1] 
b04 <- para04[2]


ntrial <- 1000     #Set the number of simulations
t_alpha <- 1/2
t_beta  <- 1/2



timestart <- Sys.time()
pp <- ECP_multiple_binary_modle(a01, b01, a02, b02, a03, b03, a04, b04, n01, x01, n02, x02, n03, x03, n04, x04, mu01, sig01, mu02, sig02, mu03, sig03, mu04, sig04, nc, pc, nt, pt, t_alpha, t_beta, ntrial)  #probability that treatment is superior to control
timeend <- Sys.time()
print(timeend - timestart)

pp_ECP_multiple_binary_1000<- data.frame(pp)
write.xlsx(pp_ECP_multiple_binary_1000, file = "pp_ECP_multiple_binary_1000.xlsx")


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
pt <- c(0.4, 0.6, 0.55, 0.6, 0.6, 0.65, 0.35, 0.4)



cl<- makeCluster(10)      
registerDoParallel(cl)      

timestart <- Sys.time()

result_pp <- foreach(i = 1:length(pc),.combine = 'cbind',.packages = c("LaplacesDemon","invgamma","runjags"))  %dopar% {
  
  ##probability that treatment is superior to control
  
  pp <- ECP_multiple_binary_modle(a01, b01, a02, b02, a03, b03, a04, b04, n01, x01, n02, x02, n03, x03, n04, x04,  mu01, sig01, mu02, sig02, mu03, sig03, mu04, sig04, nc, pc=pc[i], nt, pt=pt[i], t_alpha, t_beta, ntrial)
  
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