
###################################################
################################## We check the bootstrap coverage rates using simulation
################################## This example use the Wilcoxon R-estimator and Scheme U for bootstrap weights

######You need to install "Rccp" and "RcppArmadillo" to use "RankGarch.cpp"
install.packages("Rcpp")
install.packages("RcppArmadillo")
library(Rcpp)
library(RcppArmadillo)


sourceCpp("RankGarch.cpp")


####################################################
########################## First we simulate k replications and compute R-estimators
k= 1000
int=300 #remove first int number of initial observations to generate stationary observations
n=1000
theta = matrix(NA,3,k) #matrix to store estimators of simulations
nn=int+n
eta=rep(0,nn) #error term
x = matrix(NA, k, n)
omega=6.50E-06 #set ture value of omega using sp500 from fGarch
alpha=1.77E-01
beta=7.16E-01
theta0=c(omega,alpha,beta)
a=matrix(0,3,k) #matrix to store the iteration adjustment
weight=rep(1,n)


for(j in 1:k)
{
  ####double exponential errors
  #rv=runif(nn,0,1)
  #for(i in 1:nn)
  #{
  # if(rv[i]<0.5)
  #{eta[i]=log(2*rv[i])}
  # else
  #{eta[i]=-log(2*(1-rv[i]))}  
  #}
  #eta=eta/sqrt(2)
  
  #########normal distribution errors
  #eta=rnorm(nn,0,1)
  
  #########logistic distribution errors
  eta=rlogis(nn,0,1)
  eta=eta/1.79
  
  #######t(3) distribution
  #eta=rt(nn,df=3)
  #eta=eta/sqrt(3)
  
  #simulate the errors and observations
  VVt=XX=rep(0,nn)
  VVt[1]=omega/(1-alpha-beta)
  for(i in 2:nn)
  {
    VVt[i]=omega+alpha*XX[i-1]^2+beta*VVt[i-1]
    XX[i]=sqrt(VVt[i])*eta[i]
  }
  
  x[j, ]=XX[(1+int):nn]
  
  #use simulated data to generate R-estimator
  #########################try QMLE as the initial estimator
  fgarch_fit = fGarch::garchFit(x[j, ] ~garch(1,1), data = x[j, ], include.mean = F)
  QMLE_fGarch = fgarch_fit@fit$par ####used as the initial estimator
  
  
  theta[,j] = Resti(QMLE_fGarch, x[j, ], score = "Wilcoxon")$coef
}


theta_n = matrix(NA, 3, k)
theta_n = theta ###Note that this is R-estimator of $theta_{n\varphi}$


##################function to estimate the scalar $c_\varphi$
c_vp = function(data, theta){
  Xt_square = data^2
  Xt_square_aver = mean(Xt_square)
  
  cvp = (theta[1]/Xt_square_aver + theta[2])/(1 - theta[3])
  return(cvp)
}

#### estimated c_\{varphi}
c_Wil_Log = rep(NA, k)
for(i in 1:k){
  c_Wil_Log[i] = c_vp(x[i,], theta_n[, i])
}


################# $n^{-1/2}*(theta_{n\varphi} - theta_{0\varphi})$
sim_theta = matrix(NA, 3, k)
for(i in 1:k){
  sim_theta[1, i] = sqrt(n)*(theta_n[1, i] - theta0[1]*c_Wil_Log[i])
  sim_theta[2, i] = sqrt(n)*(theta_n[2, i] - theta0[2]*c_Wil_Log[i])
  sim_theta[3, i] = sqrt(n)*(theta_n[3, i] - theta0[3])
}



########################################### Bootstrap R-estimator and get coverage rates
B = 2000 #bootstrap size
######scheme M
#sigman=sqrt(1-1/n)
#prob=rep(1,n)/n
######scheme E
#sigman=sqrt((n-1)/(n+1)) 
######scheme U
U=runif(n,0.5,1.5)
sigman=sqrt(var(U/mean(U)))

###coverage rates
cover1 = cover2 = matrix(0, 3, k)
reject_pr = c(0.05, 0.1)

##################################################bootstrap for each replication
for(i in 1:k){
  theta_Boot = matrix(0,3,B)
  aa = matrix(0,3,B)
  Conver_Boot = rep(NA, B)
  
  for(b in 1:B){
    #scheme M
    #weight_boot = rmultinom(1, size = n, prob)
    #scheme E
    #weight_boot = rexp(n,1)
    #weight_boot = weight_boot/mean(weight_boot)
    #scheme U
    weight_boot=runif(n,0.5,1.5)
    weight_boot=weight_boot/mean(weight_boot)
    
    theta_Boot[, b] = theta_n[, i] #set initial values of Bootstrap estimate
    
    fitBoot = Boot_Resti(weight_boot, theta_Boot[, b], x[i, ], score = "Wilcoxon")
    theta_Boot[, b] = fitBoot$coef
  
  Boot_Normalized_theta = matrix(NA, 3, B)
  
  Boot_theta1_low = Boot_theta1_upper = Boot_theta2_low = Boot_theta2_upper = rep(NA, 3)
  for(l in 1:3){
  Boot_Normalized_theta[l, ] = 1/sigman*sqrt(n)*(theta_Boot[l,] - theta_n[l, i])
  
  Boot_theta1_low[l] = quantile(Boot_Normalized_theta[l, ], reject_pr[1]/2, na.rm = T)
  Boot_theta1_upper[l] = quantile(Boot_Normalized_theta[l, ], (1 - reject_pr[1]/2), na.rm = T)
  
  Boot_theta2_low[l] = quantile(Boot_Normalized_theta[l, ], reject_pr[2]/2, na.rm = T)
  Boot_theta2_upper[l] = quantile(Boot_Normalized_theta[l, ], (1 - reject_pr[2]/2), na.rm = T)
  
  if((sim_theta[l, i] >= Boot_theta1_low[l])&(sim_theta[l, i] <= Boot_theta1_upper[l])){
    cover1[l, i] = 1
  } else{cover1[l, i] = 0}
  if((sim_theta[l, i] >= Boot_theta2_low[l])&(sim_theta[l, i] <= Boot_theta2_upper[l])){
    cover2[l, i] = 1
  } else{cover2[l, i] = 0}
  }
  }
}


########coverage rates
coverage1 = coverage2 = rep(NA, 3)
for(i in 1:3){
  coverage1[i] = sum(cover1[i, ])/k
  coverage2[i] = sum(cover2[i, ])/k
}



