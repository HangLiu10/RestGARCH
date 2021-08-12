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

theta_QMLE = QMLEesti(QMLE_fGarch, r1)$coef
theta_LAD = LADesti(QMLE_fGarch, r1)$coef


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
c_LAD = c_vp(r1, theta_LAD)

vdW = theta_vdW/c(c_vdW, c_vdW, 1)
Wil = theta_Wil/c(c_Wil, c_Wil, 1)
sign_hat = theta_sign/c(c_Sign, c_Sign, 1)

QMLE = theta_QMLE/c(1, 1, 1)

LAD = theta_LAD/c(c_LAD, c_LAD, 1)


#################################
#################################
##################QQ-plot of residuals (EFCX data)
vt_vdW = vt(theta_vdW, r1)

epsilon_vdW = r1/sqrt(vt_vdW)

x= qt(ppoints(n1), df = 4.01)

par(mai = c(0.8, 0.8, 0.8, 0.8), cex.axis=1.5)

#############against t(4.01)
qqplot(x, epsilon_vdW, xlab = "t(4.01) distribution", ylab = "Residual", main = "QQ plot against t(4.01) distribution (EFCX data)")
qqline(epsilon_vdW, col = "red")

#############against t(3.01)

x= qt(ppoints(n1), df = 3.01)

par(mai = c(0.8, 0.8, 0.8, 0.8), cex.axis=1.5)

qqplot(x, epsilon_vdW, xlab = "t(3.01) distribution", ylab = "Residual", main = "QQ plot against t(3.01) distribution (EFCX data)")
qqline(epsilon_vdW, col = "red")


##################QQ-plot of residuals (S&P 500 data)
#vt_vdW = vt(theta_vdW, r1)

#epsilon_vdW = r1/sqrt(vt_vdW)

#x= qt(ppoints(n1), df = 4.01)

#par(mai = c(0.8, 0.8, 0.8, 0.8), cex.axis=1.5)

#############against t(4.01)
#qqplot(x, epsilon_vdW, xlab = "t(4.01) distribution", ylab = "Residual", main = "QQ plot against t(4.01) (S&P 500 data)")
#qqline(epsilon_vdW, col = "red")

#############against t(6.01)

#x= qt(ppoints(n1), df = 6.01)

#par(mai = c(0.8, 0.8, 0.8, 0.8), cex.axis=1.5)

#qqplot(x, epsilon_vdW, xlab = "t(6.01) distribution", ylab = "Residual", main = "QQ plot against t(3.01) (S&P 500 data)")
#qqline(epsilon_vdW, col = "red")







