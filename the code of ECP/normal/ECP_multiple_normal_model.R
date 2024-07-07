###################################ECP_multiple_normal########################################
##load packages
library(dplyr)
library(LaplacesDemon)
library(invgamma)
library("foreach")
library("doParallel")
library(openxlsx)


#-----Function to decide tuning parameter a and b in elastic function g(T)=exp⁡(a+b∙log⁡(T))------------#
# Inputs:
# x0: external data
# n0: sample size for external data
# nc: sample size for control arms
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


#Defeine JAGS model
ECP_normal_multiple_modle =  "
model {

  #current likelihood
  for (i in 1:nc) {
    xc[i] ~ dnorm(muc, sig)
  }
  
  muc ~ dnorm(mu_sum, sig_sum)
  sig ~ dgamma(0.001,0.001)
} 
"
writeLines(ECP_normal_multiple_modle,con = "ECP_normal_multiple_modle.txt")


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

ECP_multiple_normal_model <- function(a01, b01, a02, b02, a03, b03, a04, b04, n01, x01, n02, x02, n03, x03, n04, x04, nc, uc, sigc, nt, ut, sigt, ntrial){
  pp <- numeric(ntrial)
  mu01 <- mean(x01)
  mu02 <- mean(x02)
  mu03 <- mean(x03)
  mu04 <- mean(x04)
  
  sig01 <- var(x01)
  sig02 <- var(x02)
  sig03 <- var(x03)
  sig04 <- var(x04)
  
  for(trial in 1:ntrial){
    set.seed(100 + trial)
    # generate control and treatment data
    xc <- rnorm(nc, uc, sigc)
    xt <- rnorm(nt, ut, sigt)
    # calculate statistic and g(t)
    sp01 <- ((n01-1)*var(x01) + (nc-1)*var(xc))/(n01 + nc - 2) # pooled variance
    T01 <- max(n01, nc)^(-1/4)*abs(mean(x01)-mean(xc))/(sqrt(sp01/n01 + sp01/nc))
    gt01 <- exp(a01 + b01*log(T01))
    if(gt01 < 0.00001){gt01 <- 0.00001}
    
    sp02 <- ((n02-1)*var(x02) + (nc-1)*var(xc))/(n02 + nc - 2) # pooled variance
    T02 <- max(n02, nc)^(-1/4)*abs(mean(x02)-mean(xc))/(sqrt(sp02/n02 + sp02/nc))
    gt02 <- exp(a02 + b02*log(T02))
    if(gt02 < 0.00001){gt02 <- 0.00001}
    
    sp03 <- ((n03-1)*var(x03) + (nc-1)*var(xc))/(n03 + nc - 2) # pooled variance
    T03 <- max(n03, nc)^(-1/4)*abs(mean(x03)-mean(xc))/(sqrt(sp03/n03 + sp03/nc))
    gt03 <- exp(a03 + b03*log(T03))
    if(gt03 < 0.00001){gt03 <- 0.00001}
    
    sp04 <- ((n04-1)*var(x04) + (nc-1)*var(xc))/(n04 + nc - 2) # pooled variance
    T04 <- max(n04, nc)^(-1/4)*abs(mean(x04)-mean(xc))/(sqrt(sp04/n04 + sp04/nc))
    gt04 <- exp(a04 + b04*log(T04))
    if(gt04 < 0.00001){gt04 <- 0.00001}
    
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
    monitor = c("muc","sig")
    
    ##Run JAGS
    results <- run.jags(model = "ECP_normal_multiple_modle.txt",
                        monitor = monitor,
                        data = dat,
                        n.chains = 3,
                        inits = c(inits1,inits2,inits3),
                        plots = FALSE,
                        burnin = 10000,
                        sample = 10000,
                        thin = 1)
    chains = rbind(results$mcmc[[1]],results$mcmc[[2]],results$mcmc[[3]])
    muc = sample(chains[,"muc"],10000)
    
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
## Generate four simulated external data
set.seed(1)
u0 <- 1
sig <- 1
mu <- round(rnorm(4, u0, 0.1), 2)

set.seed(14767)
# generate external data1
n01 <- 40
x01 <- rnorm(n01, mu[1], sig)

set.seed(13501)
# generate external data2
n02 <- 50
x02 <- rnorm(n02, mu[2], sig)

set.seed(4640)
# generate external data3
n03 <- 45
x03 <- rnorm(n03, mu[3], sig)

set.seed(2968)
# generate external data4
n04 <- 55
x04 <- rnorm(n04, mu[4], sig)



## control and treatment
nc <- 25 
nt <- 50

uc <- (mean(x01)+mean(x02)+mean(x03)+mean(x04))/4
ut <- (mean(x01)+mean(x02)+mean(x03)+mean(x04))/4

sigc <- sqrt((var(x01)*(n01-1) + var(x02)*(n02-1) + var(x03)*(n03-1) + var(x04)*(n04-1))/(n01+n02+n03+n04-4))
sigt <- sqrt((var(x01)*(n01-1) + var(x02)*(n02-1) + var(x03)*(n03-1) + var(x04)*(n04-1))/(n01+n02+n03+n04-4))



#Determine q1 and q2 based on the grid search results
para01 <- decide_para(x0 = x01, n0 = n01, nc, gamma=1, q1=0.93, q2=0.01, small = 0.00000000001, large = 100000000000, R = 50000)
para02 <- decide_para(x0 = x02, n0 = n02, nc, gamma=1, q1=0.93, q2=0.01, small = 0.00000000001, large = 100000000000, R = 50000)
para03 <- decide_para(x0 = x03, n0 = n03, nc, gamma=1, q1=0.93, q2=0.01, small = 0.00000000001, large = 100000000000, R = 50000)
para04 <- decide_para(x0 = x04, n0 = n04, nc, gamma=1, q1=0.93, q2=0.01, small = 0.00000000001, large = 100000000000, R = 50000)

a01 <- para01[1] 
b01 <- para01[2]
a02 <- para02[1] 
b02 <- para02[2]
a03 <- para03[1] 
b03 <- para03[2]
a04 <- para04[1] 
b04 <- para04[2]


ntrial <- 1000     #Set the number of simulations
timestart <- Sys.time()
pp <- ECP_multiple_normal_model(a01, b01, a02, b02, a03, b03, a04, b04, n01, x01, n02, x02, n03, x03, n04, x04, nc, uc, sigc, nt, ut, sigt, ntrial)
timeend <- Sys.time()

print(timeend - timestart)

pp_ECP_multiple_normal_1000<- data.frame(pp)
write.xlsx(pp_ECP_multiple_normal_1000, file = "pp_ECP_multiple_normal_1000.xlsx")


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
  
  pp <- ECP_multiple_normal_model(a01, b01, a02, b02, a03, b03, a04, b04, n01, x01, n02, x02, n03, x03, n04, x04, nc, uc=uc[i], sigc, nt, ut=ut[i], sigt, ntrial)
  
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