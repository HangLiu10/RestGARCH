
##################################################
######################################Simulation study of the R-estimators and QMLE
######################################Data generating process: GARCH(2, 1)

######You need to install "Rccp" and "RcppArmadillo" to use "RankGarch.cpp"
install.packages("Rcpp")
install.packages("RcppArmadillo")
library(Rcpp)
library(RcppArmadillo)



sourceCpp("RankGarch.cpp")



#############Simluate k samples for GARCH(2, 1) with various error distributions
k=500 #num of replications
int=300 #remove first int number of initial observations to generate stationary observations
n=1000 #sample size
theta_sign_sim = theta_vdW_sim = theta_Wil_sim = theta_QMLE_sim = matrix(NA,4,k) #matrix to store estimators of simulations
nn=int+n
eta=rep(0,nn) #error term
x = matrix(NA, k, n)
omega= 4.46E-06
alpha1= 0.0525
alpha2= 0.108
beta= 0.832
theta0=c(omega,alpha1,alpha2, beta)
sign_change = vdW_change = Wil_change = QMLE_change = matrix(0,4,k) #matrix to store the iteration adjustment



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
  eta=rnorm(nn,0,1)
  
  #########logistic distribution errors
  #eta=rlogis(nn,0,1)
  #eta=eta/1.79
  
  #######t(3) distribution
  #eta=rt(nn,df=3)
  #eta=eta/sqrt(3)
  
  #simulate the errors and observations
  VVt=XX=rep(0,nn)
  VVt[1]=VVt[2]=omega/(1-alpha1-alpha2-beta)
  for(i in 3:nn)
  {
    VVt[i]=omega+alpha1*XX[i-1]^2+alpha2*XX[i-2]^2+beta*VVt[i-1]
    XX[i]=sqrt(VVt[i])*eta[i]
  }
  
  x[j, ]=XX[(1+int):nn]
}



#######################R-estimation and QMLE for the simulated samples
for(j in 1:k){
  ############################use simulated data to generate R-estimator
  #########################try QMLE as the initial estimator
  fgarch_fit = fGarch::garchFit(x[j, ] ~garch(2,1), data = x[j, ], include.mean = F)
  QMLE_fGarch = fgarch_fit@fit$par ####used as the initial estimator
  
  theta_sign = QMLE_fGarch #the initial values in iterations (one may try different values)
  theta_vdW = QMLE_fGarch
  theta_Wil = QMLE_fGarch
  theta_QMLE = QMLE_fGarch
  
  #theta_sign = theta_Wil = theta_vdW = theta0
  #theta_QMLE = theta0
  
  #########sign
  fit_sign = Resti_21(theta_sign, x[j, ], score = "Sign")
  sign_change[, j] = fit_sign$change
  theta_sign = fit_sign$coef
  
  #########vdW
  fit_vdW = Resti_21(theta_vdW, x[j, ], score = "vdW")
  vdW_change[, j] = fit_vdW$change
  theta_vdW = fit_vdW$coef
  
  #########Wilcoxon
  fit_Wil = Resti_21(theta_Wil, x[j, ], score = "Wilcoxon")
  Wil_change[, j] = fit_Wil$change
  theta_Wil = fit_Wil$coef
  
  #############QMLE
  fit_QMLE = QMLEesti_21(theta_QMLE, x[j, ])
  QMLE_change[, j] = fit_QMLE$change
  theta_QMLE = fit_QMLE$coef
  
  
  
  theta_sign_sim[,j] = theta_sign
  theta_vdW_sim[,j] = theta_vdW
  theta_Wil_sim[,j] = theta_Wil
  theta_QMLE_sim[,j] = theta_QMLE
}


################################################compute bias and MSE

##################function to estimate the scalar $c_\varphi$
c_vp = function(data, theta){
  Xt_square = data^2
  Xt_square_aver = mean(Xt_square)
  
  cvp = (theta[1]/Xt_square_aver + theta[2] + theta[3])/(1 - theta[4])
  return(cvp)
}



################################################compute bias and MSE
#### estimated c_\{varphi}
c_vdW = c_Wil = c_Sign = rep(NA, k)
for(i in 1:k){
  c_vdW[i] = c_vp(x[i,], theta_vdW_sim[, i])
  c_Wil[i] = c_vp(x[i,], theta_Wil_sim[, i])
  c_Sign[i] = c_vp(x[i,], theta_sign_sim[, i])
}


##################bias and MSE
bias_Wil = bias_Sign = bias_vdW = bias_QMLE = MSE_Wil = MSE_Sign = MSE_vdW = MSE_QMLE = rep(NA, 4)
c1 = c2 = c3 = matrix(NA, 4, k)
for(j in 1:k){
  c1[, j] = c(c_Wil[j], c_Wil[j], c_Wil[j], 1)
  c2[, j] = c(c_Sign[j], c_Sign[j], c_Sign[j], 1)
  c3[, j] = c(c_vdW[j], c_vdW[j], c_vdW[j], 1)
}


for(i in 1:4){
  bias_Wil[i] = mean(theta_Wil_sim[i, ]/c1[i] - theta0[i], na.rm = T)
  MSE_Wil[i] = mean((theta_Wil_sim[i, ]/c1[i] - theta0[i])^2, na.rm = T)
  
  bias_Sign[i] = mean(theta_sign_sim[i, ]/c2[i] - theta0[i], na.rm = T)
  MSE_Sign[i] = mean((theta_sign_sim[i, ]/c2[i] - theta0[i])^2, na.rm = T)
  
  bias_vdW[i] = mean(theta_vdW_sim[i, ]/c3[i] - theta0[i], na.rm = T)
  MSE_vdW[i] = mean((theta_vdW_sim[i, ]/c3[i] - theta0[i])^2, na.rm = T)
  
  bias_QMLE[i] = mean(theta_QMLE_sim[i, ] - theta0[i], na.rm = T)
  MSE_QMLE[i] = mean((theta_QMLE_sim[i, ] - theta0[i])^2, na.rm = T)
}



