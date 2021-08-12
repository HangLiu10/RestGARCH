##################################################Boxplots to produce Figure 1
######################################Simulation study of the vdW R-estimator and QMLE
######################################Data generating process: GARCH(1, 1)

######You need to install "Rcpp" and "RcppArmadillo" to use "RankGarch.cpp"
install.packages("Rcpp")
install.packages("RcppArmadillo")
library(Rcpp)
library(RcppArmadillo)



sourceCpp("RankGarch.cpp")

install.packages("fGarch") #get QMLE from "fGarch"
library(fGarch)

install.packages("sn") # generate skew-t dist.
library(sn)


#############Simluate k samples for GARCH(1, 1) with various error distributions
k=1000 #num of replications
int=300 #remove first int number of initial observations to generate stationary observations
n=1000 #sample size
theta_vdW_sim = theta_QMLE_sim = matrix(NA,3,k) #matrix to store estimators of simulations
nn=int+n
eta=rep(0,nn) #error term
Xt = matrix(NA, k, n)  #observations
omega=6.50E-06 #set ture value of omega using sp500 data from fGarch
alpha=1.77E-01
beta=7.16E-01
theta0=c(omega,alpha,beta)
vdW_change = QMLE_change = matrix(0,3,k) #matrix to store the iteration adjustment

for(j in 1:k){
  set.seed(j)
  
  #########normal distribution errors
  eta=rnorm(nn,0,1)
  
  ##########skew normal
  #########################generate a large sample to compute the mean of skew-normal dist.
  #mean_skewnorm = mean(rst(10000, xi= 0, omega = (13*pi)/(13*pi - 25), alpha = 5, nu=Inf))
  #eta = rst(nn, xi= 0, omega = (13*pi)/(13*pi - 25), alpha = 5, nu=Inf)
  #eta = (eta - mean_skewnorm)/sd(eta)
 
  #simulate the errors and observations
  VVt=XX=rep(0,nn)
  VVt[1]=omega/(1-alpha-beta)
  for(i in 2:nn)
  {
    VVt[i]=omega+alpha*XX[i-1]^2+beta*VVt[i-1]
    XX[i]=sqrt(VVt[i])*eta[i]
  }
  
  Xt[j, ]=XX[(1+int):nn]
}

#######################R-estimators (of $\theta_{0\varphi}$) and QMLE for the simulated samples
for(j in 1:k){
  ############################use simulated data to generate R-estimator
  #########################try QMLE as the initial estimator
  fgarch_fit = fGarch::garchFit(Xt[j, ] ~garch(1,1), data = Xt[j, ], include.mean = F)
  QMLE_fGarch = fgarch_fit@fit$par ####used as the initial estimator
  
  theta_vdW = QMLE_fGarch #the initial values in iterations (one may try different values)
  theta_QMLE = QMLE_fGarch
  
  #theta_vdW = theta0
  #theta_QMLE = theta0
  
  #########vdW
  fit_vdW = Resti(theta_vdW, Xt[j, ], score = "vdW")
  vdW_change[, j] = fit_vdW$change
  theta_vdW = fit_vdW$coef

  ##############QMLE
  fit_QMLE = QMLEesti(theta_QMLE, Xt[j, ])
  QMLE_change[, j] = fit_QMLE$change
  theta_QMLE = fit_QMLE$coef

  theta_vdW_sim[,j] = theta_vdW
  theta_QMLE_sim[,j] = theta_QMLE
}

##################function to estimate the scalar $c_\varphi$
c_vp = function(data, theta){
  Xt_square = data^2
  Xt_square_aver = mean(Xt_square)
  
  cvp = (theta[1]/Xt_square_aver + theta[2])/(1 - theta[3])
  return(cvp)
}



################################################compute bias and MSE
#### estimated c_\{varphi}
c_vdW = rep(NA, k)
for(i in 1:k){
  c_vdW[i] = c_vp(Xt[i,], theta_vdW_sim[, i])
}




#######################################bias and MSE
bias_vdW = bias_QMLE = MSE_vdW = MSE_QMLE = rep(NA, 3)
c3 = matrix(NA, 3, k)
for(j in 1:k){
  c3[, j] = c(c_vdW[j], c_vdW[j], 1)
}


###########R-estimator of true parameter
vdW_hat = matrix(NA,3,k)

for(i in 1:3){
  vdW_hat[i, ] = theta_vdW_sim[i, ]/c3[i, ]
}


for(i in 1:3){
  bias_vdW[i] = mean(theta_vdW_sim[i, ]/c3[i, ] - theta0[i], na.rm = T)
  MSE_vdW[i] = mean((theta_vdW_sim[i, ]/c3[i, ] - theta0[i])^2, na.rm = T)
  
  bias_QMLE[i] = mean(theta_QMLE_sim[i, ] - theta0[i], na.rm = T)
  MSE_QMLE[i] = mean((theta_QMLE_sim[i, ] - theta0[i])^2, na.rm = T)
}


##################We will also include boxplot of R-estimator of Andrews (2012)
#########This part of codes illustrate how to obtain R-estimator of Andrews (2012)
#-------------------------------------------------
# Choose GARCH model orders
#--------------------------------------------------
p<- 1
q<- 1


###################################################

#---------------------------------------------------------------------------------
#Auxiliary functions:
#---------------------------------------------------------------------------------

###################################################
#--------------------------------------------------
# the weight function lambda.
# Here we have a weight function which is optimal when the GARCH noise is rescaled student's t with df=7
#-------------------------------------------------
#lambda <- function(x){
#  if(x<1) {c<-(7*(sqrt(5/7)*qt((x+1)/2,7))^2-5)/((sqrt(5/7)*qt((x+1)/2,7))^2+5)}
#  if(x==1) {c<-7}
#  return(c)
#}
###################################################

#vdW score lambda
lambda <- function(x){
  c = (qnorm((x+1)/2))^2 -1
  return(c)
}

###################################################
#--------------------------------------------------
# the weight function lambda squared.
#-------------------------------------------------
lambda2 <- function(x){
  c <- lambda(x)^2
  return(c)
}
###################################################
theta_hat_BA = matrix(NA,2,k)

for(l in 1:k){
  X = Xt[l, ]
  ###################################################
  #--------------------------------------------------
  # getD returns the value of D for a candidate parameter vector theta
  #-------------------------------------------------
  getD <- function(theta){
    #--------------------------------------------------
    # calculate sigma square
    sigmasq <- array(1, (n+max(0,q-p)))
    for (t in (p+1):n){
      sigmasq[t+max(0,q-p)]<-1+sum(theta[1:p]*((X[(t-1):(t-p)])^2))
      if(q>0){sigmasq[t+max(0,q-p)]<-sigmasq[t+max(0,q-p)]+sum(theta[(p+1):(p+q)]*sigmasq[(t+max(0,q-p)-1):(t+max(0,q-p)-q)])}
    }
    
    #--------------------------------------------------
    # calculate epsilon
    epsilon <- log(X^2)-log(abs(sigmasq[(1+max(0,q-p)):(n+max(0,q-p))]))
    #--------------------------------------------------
    # sort epsilon
    sortep <- sort(epsilon[p+1:n])
    #--------------------------------------------------
    # calculate D
    D<-0
    for (t in 1:(n-p)){
      D<-D+lambda(t/(n-p+1))*(sortep[t]-mean(sortep))
    }
    if(min(theta)<0){D=D+10000000}
    if((q>0) & (sum(theta[(p+1):(p+q)])>=1)){D=D+10000000}
    return(D)
  }
  
  
  #---------------------------------------------------------------------------------
  #End of auxiliary functions
  #---------------------------------------------------------------------------------
  
  ###################################################
  #--------------------------------------------------
  # Main Body
  
  #--------------------------------------------------
  # generate 100 initial values for theta=(alpha1/alpha0, ... , alphap/alpha0, beta1, ... ,betaq)
  # under the condition that the sum of alpha1, ... ,alphap, beta1, ... ,betaq is less than 1
  #--------------------------------------------------
  no_initial <- 100
  temp_alphabeta <-matrix(0,nr=no_initial,nc=p+q+1)
  temp_theta <-matrix(0,nr=no_initial,nc=p+q)
  count<-0
  repeat{
    alphabeta <- runif(p+q)
    if (sum(alphabeta)<1){
      count<-count + 1
      temp_alphabeta[count,] <- c(var(X)*(1-sum(alphabeta)),alphabeta)
    }
    if(count == no_initial)break 
  }
  
  for (i in 1:p){
    temp_theta[,i]<-temp_alphabeta[,(i+1)]/temp_alphabeta[,1]
  }
  if (q>0){
    for (i in (p+1):(p+q)){
      temp_theta[,i]<-temp_alphabeta[,(i+1)]
    }
  }
  
  #--------------------------------------------------
  # calculate D for the initial parameter vectors theta
  #--------------------------------------------------
  temp_D<- array(1000,no_initial)
  for (i in 1:no_initial){
    temp_D[i] <-getD(temp_theta[i,])
  }
  
  #--------------------------------------------------
  # find the 3 initial values with the smallest values for D
  #--------------------------------------------------
  store_D <- array(1000,3)
  store_theta <-matrix(0,nr=3,nc=p+q)
  for (j in 1:3){
    store_D[j] <- min(temp_D)
    store_theta[j,] <- temp_theta[which.min(temp_D),]
    temp_D[which.min(temp_D)] <- max(temp_D)+1
  }
  
  #--------------------------------------------------
  # Using the 3 initial values as starting points, find optimized values for theta
  # The Rank estimate has the smallest corresponding D value
  #--------------------------------------------------
  min_theta <- store_theta[1,]
  min_D <- store_D[1]
  for (j in 1:3){
    theta<-store_theta[j,]
    theta<-optim(par=theta,fn=getD)
    tempD<-theta$value
    if (tempD<min_D){
      min_D<-tempD
      min_theta<-theta$par
    }
  }
  
  #--------------------------------------------------
  # Give Rank estimate of theta
  #--------------------------------------------------
  cat ("theta=(alpha1/alpha0,...,alphap/alpha0,beta1,...,betaq):",min_theta,"\n")
  theta_hat_BA[,l] = min_theta
  print(l)
}
#--------------------------------------------------


#save(theta_hat_BA, file = "AndrewsNormal.RData")




########## estimate $\omega$ and compute the bias and MSE of the estimator

##################function to estimate $\omega$
c_BA = function(data, theta_BA){
  Xt_square = data^2
  Xt_square_aver = mean(Xt_square)
  
  cBA = ((1 - theta_BA[2])*Xt_square_aver)/(1 + theta_BA[1]*Xt_square_aver)
  return(cBA)
}


omega_BA = rep(NA, k)

for(i in 1:k){
  omega_BA[i] = c_BA(data = Xt[i,], theta_BA = theta_hat_BA[,i])
}

hat_BA = matrix(NA, 3, k)

hat_BA[1, ] = omega_BA
hat_BA[2, ] = omega_BA*theta_hat_BA[1,]
hat_BA[3, ] = theta_hat_BA[2,]


bias_BA = MSE_BA = rep(NA, 3)


for(i in 1:3){
  bias_BA[i] = mean(hat_BA[i, ] - theta0[i], na.rm = T)
  MSE_BA[i] = mean((hat_BA[i, ] - theta0[i])^2, na.rm = T)
}






###########################boxplot of the vdW, QMLE and Andrews' R-estimator
vdW = vdW_hat
QMLE = theta_QMLE_sim
BA = hat_BA

R = k
#####get values for boxplot
names=c(rep("QMLE", R), rep("vdW(BA)", R), rep("vdW", R))
value=c(QMLE[1,],  BA[1, ], vdW[1, ])
data=data.frame(names,value)

data$names = factor(data$names, levels = c("QMLE",  "vdW(BA)", "vdW"))


#Draw the boxplot, with the number of individuals per group
par(mfcol = c(1,3), mar = c(2.2, 2.2, 2.2, 2.2), cex.axis = 1.5)

######boxplot for omega
a=boxplot(data$value~data$names, main = expression(omega), lwd = 0.8, lex.order = F, cex.main = 2, ylim = c(0, 2.2E-05))
abline(h = omega, col = "red")
MSEratio1 = MSE_QMLE[1]/MSE_BA[1]
MSEratio1 = format(round(MSEratio1, 3), nsmall = 3)

MSEratio2 = MSE_QMLE[1]/MSE_vdW[1]
MSEratio2 = format(round(MSEratio2, 3), nsmall = 3)

text(x= 1, y= 1.8E-05, labels= "")
text(x= 2, y= 1.8E-05, labels= MSEratio1, cex = 1.5)
text(x= 3, y= 1.8E-05, labels= MSEratio2, cex = 1.5)

#####boxplot for alpha
value=c(QMLE[2,],  BA[2, ], vdW[2, ])

data=data.frame(names,value)

data$names = factor(data$names, levels = c("QMLE",  "vdW(BA)", "vdW"))


a=boxplot(data$value~data$names, main = expression(alpha), lwd = 0.8, lex.order = F, cex.main = 2, ylim = c(0.06, 0.34))
abline(h = alpha, col = "red")
MSEratio1 = MSE_QMLE[2]/MSE_BA[2]
MSEratio1 = format(round(MSEratio1, 3), nsmall = 3)

MSEratio2 = MSE_QMLE[2]/MSE_vdW[2]
MSEratio2 = format(round(MSEratio2, 3), nsmall = 3)


text(x= 1, y= 0.29, labels= "")
text(x= 2, y= 0.29, labels= MSEratio1, cex = 1.5)
text(x= 3, y= 0.29, labels= MSEratio2, cex = 1.5)


#############boxplot for beta

value=c(QMLE[3,],  BA[3, ], vdW[3, ])

data=data.frame(names,value)

data$names = factor(data$names, levels = c("QMLE",  "vdW(BA)", "vdW"))


a=boxplot(data$value~data$names, main = expression(beta), lwd = 0.8, lex.order = F, cex.main = 2, ylim = c(0.42, 0.89))
abline(h = beta, col = "red")
MSEratio1 = MSE_QMLE[3]/MSE_BA[3]
MSEratio1 = format(round(MSEratio1, 3), nsmall = 3)

MSEratio2 = MSE_QMLE[3]/MSE_vdW[3]
MSEratio2 = format(round(MSEratio2, 3), nsmall = 3)


text(x= 1, y= 0.8, labels= "")
text(x= 2, y= 0.8, labels= MSEratio1, cex = 1.5)
text(x= 3, y= 0.8, labels= MSEratio2, cex = 1.5)



