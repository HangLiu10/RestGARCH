######You need to install "Rccp" and "RcppArmadillo" to use "RankGarch.cpp"
install.packages("Rcpp")
install.packages("RcppArmadillo")
library(Rcpp)
library(RcppArmadillo)



sourceCpp("RankGarch.cpp")


####################read data
data1=read.csv("EFCX.csv",header=T,sep=",",stringsAsFactors = F)
price=data1$Adj.Close  #daily price of EFCX
#price = data1$GBP.USD
len1 = length(price)
date1=as.Date(data1$Date[1:len1],"%d/%m/%Y")
plot(date1, price, type = "l", main="EFCX Daily Price",xlab = "Date",ylab="Daily Price")



#log-returns
r1=log(price[2:len1]/price[1:(len1-1)])
n1=length(r1)
date1=as.Date(data1$Date[2:len1],"%d/%m/%Y")

plot(date1,r1,type = "l",main="Returns of EFCX",xlab = "Date",ylab="Return")


########estimate by fGarch
library(fGarch)
fit1=garchFit(r1~garch(1,1),data=r1,trace=F,include.mean = F)
QMLE_fGarch = fit1@fit$par

#########R-estimators and QMLE
theta_vdW = Resti(QMLE_fGarch, r1, score = "vdW")$coef
theta_sign = Resti(QMLE_fGarch, r1, score = "Sign")$coef
theta_Wil = Resti(QMLE_fGarch, r1, score = "Wilcoxon")$coef


##################function to estimate $c_\varphi$
c_vp = function(data, theta){
  Xt_square = data^2
  Xt_square_aver = mean(Xt_square)
  
  cvp = (theta[1]/Xt_square_aver + theta[2])/(1 - theta[3])
  return(cvp)
}

###compare with the estimated c_vp
c_vdW = c_vp(r1, theta_vdW)
c_Wil = c_vp(r1, theta_Wil)
c_Sign = c_vp(r1, theta_sign)

vdW = theta_vdW/c(c_vdW, c_vdW, 1)
Wil = theta_Wil/c(c_Wil, c_Wil, 1)
sign_hat = theta_sign/c(c_Sign, c_Sign, 1)

