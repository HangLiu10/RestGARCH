#/*#include <Rcpp.h>*/
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
NumericVector vt(NumericVector theta, NumericVector r){
  int n = r.size();
  NumericVector v(n);

  v[0] = theta[0]/(1-theta[2]);
  for(int i = 1; i < n; ++i){
    v[i] = theta[0] + theta[1]*pow(r[i-1],2) + theta[2]*v[i-1];
  }

  return(v);
}


// [[Rcpp::export]]
NumericVector vt_21(NumericVector theta, NumericVector r){ /*vt for GARCH(2, 1)*/
  int n = r.size();
  NumericVector v(n);

  v[0] = theta[0]/(1-theta[3]);
  v[1] = theta[0] + theta[1]*pow(r[0], 2) + theta[3]*v[0];

  for(int i = 2; i < n; ++i){
    v[i] = theta[0] + theta[1]*pow(r[i-1], 2) + theta[2]*pow(r[i-2], 2) + theta[3]*v[i-1];
  }

  return(v);
}



NumericMatrix der(NumericVector theta, NumericVector r){
  int n = r.size();
  NumericMatrix x(3, n);
  double b = 1 - theta[2];
  NumericVector v = vt(theta, r);

  x(0,0) = 1/b;
  x(1,0) = 0;
  x(2,0) = theta[0]/pow(b, 2);

  for(int i = 1; i < n; ++i){
    x(0, i) = 1 + theta[2]*x(0, i - 1);
    x(1, i) = pow(r[i - 1], 2) + theta[2]*x(1, i - 1);
    x(2, i) = v[i - 1] + theta[2]*x(2, i - 1);
  }

  return(x);
}



NumericMatrix der_21(NumericVector theta, NumericVector r){ /* \dot{v}t for GARCH(2, 1)*/
  int n = r.size();
  NumericMatrix x(4, n);
  double b = 1 - theta[3];
  NumericVector v = vt_21(theta, r);

  x(0,0) = 1/b;
  x(1,0) = 0;
  x(2,0) = 0;
  x(3,0) = theta[0]/pow(b, 2);

  x(0, 1) = 1 + theta[3]*x(0, 0);
  x(1, 1) = pow(r[0], 2) + theta[3]*x(1, 0);
  x(2, 1) = 0;
  x(3, 1) = v[0] + theta[3]*x(3, 0);

  for(int i = 2; i < n; ++i){
    x(0, i) = 1 + theta[3]*x(0, i - 1);
    x(1, i) = pow(r[i - 1], 2) + theta[3]*x(1, i - 1);
    x(2, i) = pow(r[i - 2], 2) + theta[3]*x(2, i - 1);
    x(3, i) = v[i - 1] + theta[3]*x(3, i - 1);
  }

  return(x);
}




/*R-estimators for GARCH(1,1) and GARCH(2,1); compute change at each step*/
NumericVector R_change(NumericVector theta, NumericVector r, String score){
  int n = r.size();
  NumericVector v = vt(theta, r);
  NumericVector ran(n), w(n);
  NumericMatrix x(3, n);

  x = der(theta,r);
  w = 1/pow(v, 2);

  NumericVector err = r/pow(v, 0.5);
  NumericVector sorted = clone(err).sort();
  ran = match(err, sorted);
  NumericVector uu = ran/(n + 1);
  NumericVector u = uu - 0.5;
  NumericVector sgn(n);
  NumericVector sgn1(n);

  sgn = ifelse(u > 0, 1, -1);
  ran = sgn*(score=="Sign") + u*(score=="Wilcoxon") + qnorm(uu)*(score=="vdW");
  /*ran = qnorm(uu)*(score=="vdW") ;*/

  arma::mat a(3, 3), b(3, 1); /*define two matrices*/
  a.zeros();
  b.zeros();

  arma::mat xx = as<arma::mat>(x); /*covert to arma mat*/
  arma::vec xcol(3);

  for(int i = 0; i < n; ++i){
    xcol = xx.col(i);
       a = a + w[i]*(xcol*xcol.t());
       for(int j = 0; j < 3; ++j){
           b(j, 0) = b(j, 0) + xx(j, i)/v[i]*(1 - ran[i]*err[i]);
           /*b(j, 0) = b(j, 0) + xx(j, i)/v[i]*(1 - ran[i]*ran[i]);*/
       }
  }

  a = inv(a);
  return(wrap(a*b));
}



NumericVector R_change_21(NumericVector theta, NumericVector r, String score){
  int n = r.size();
  NumericVector v = vt_21(theta, r);
  NumericVector ran(n), w(n);
  NumericMatrix x(4, n);

  x = der_21(theta,r);
  w = 1/pow(v, 2);

  NumericVector err = r/pow(v, 0.5);
  NumericVector sorted = clone(err).sort();
  ran = match(err, sorted);
  NumericVector uu = ran/(n + 1);
  NumericVector u = uu - 0.5;
  NumericVector sgn(n);

  sgn = ifelse(u > 0, 1, -1);
  ran = sgn*(score=="Sign") + u*(score=="Wilcoxon") + qnorm(uu)*(score=="vdW");

  arma::mat a(4, 4), b(4, 1); /*define two matrices*/
  a.zeros();
  b.zeros();

  arma::mat xx = as<arma::mat>(x); /*covert to arma mat*/
  arma::vec xcol(4);

  for(int i = 0; i < n; ++i){
    xcol = xx.col(i);
       a = a + w[i]*(xcol*xcol.t());
       for(int j = 0; j < 4; ++j){
           b(j, 0) = b(j, 0) + xx(j, i)/v[i]*(1 - ran[i]*err[i]);
       }
  }

  a = inv(a);
  return(wrap(a*b));
}



/*R-estimators for GARCH(1,1)*/
// [[Rcpp::export]]
List Resti(NumericVector theta, NumericVector r, String score, int iter = 25, double StSize = 1){
  int dim = theta.size();
  NumericVector a(dim), the(dim);

  the = clone(theta);

  for(int i = 0; i < iter; ++i){
    a = R_change(the, r, score);
    the = the - a/StSize;
  }

  return List::create(
                      _["change"] = a,
                      _["coef"] = the);
}


/*R-estimators for GARCH(2,1)*/
// [[Rcpp::export]]
List Resti_21(NumericVector theta, NumericVector r, String score, int iter = 25, double StSize = 1){
  int dim = theta.size();
  NumericVector a(dim), the(dim);

  the = clone(theta);

  for(int i = 0; i < iter; ++i){
    a = R_change_21(the, r, score);
    the = the - a/StSize;
  }

  return List::create(
                      _["change"] = a,
                      _["coef"] = the);
}




/*Bootstrap estimators for GARCH(1,1) and GARCH(2,1)*/
NumericVector Boot_R_change(NumericVector weight, NumericVector theta, NumericVector r, String score){
  int n = r.size();
  NumericVector v = vt(theta, r);
  NumericVector ran(n), w(n);
  NumericMatrix x(3, n);

  x = der(theta,r);
  w = 1/pow(v, 2);

  NumericVector err = r/pow(v, 0.5);
  NumericVector sorted = clone(err).sort();
  ran = match(err, sorted);
  NumericVector uu = ran/(n + 1);
  NumericVector u = uu - 0.5;
  NumericVector sgn(n);

  sgn = ifelse(u > 0, 1, -1);
  ran = sgn*(score=="Sign") + u*(score=="Wilcoxon") + qnorm(uu)*(score=="vdW");

  arma::mat a(3, 3), b(3, 1); /*define two matrices*/
  a.zeros();
  b.zeros();

  arma::mat xx = as<arma::mat>(x); /*covert to arma mat*/
  arma::vec xcol(3);

  for(int i = 0; i < n; ++i){
    xcol = xx.col(i);
       a = a + weight[i]*w[i]*(xcol*xcol.t());
       for(int j = 0; j < 3; ++j){
           b(j, 0) = b(j, 0) + weight[i]*xx(j, i)/v[i]*(1 - ran[i]*err[i]);
       }
  }

  a = inv(a);
  return(wrap(a*b));
}


NumericVector Boot_R_change_21(NumericVector weight, NumericVector theta, NumericVector r, String score){
  int n = r.size();
  NumericVector v = vt_21(theta, r);
  NumericVector ran(n), w(n);
  NumericMatrix x(4, n);

  x = der_21(theta,r);
  w = 1/pow(v, 2);

  NumericVector err = r/pow(v, 0.5);
  NumericVector sorted = clone(err).sort();
  ran = match(err, sorted);
  NumericVector uu = ran/(n + 1);
  NumericVector u = uu - 0.5;
  NumericVector sgn(n);

  sgn = ifelse(u > 0, 1, -1);
  ran = sgn*(score=="Sign") + u*(score=="Wilcoxon") + qnorm(uu)*(score=="vdW");

  arma::mat a(4, 4), b(4, 1); /*define two matrices*/
  a.zeros();
  b.zeros();

  arma::mat xx = as<arma::mat>(x); /*covert to arma mat*/
  arma::vec xcol(4);

  for(int i = 0; i < n; ++i){
    xcol = xx.col(i);
       a = a + weight[i]*w[i]*(xcol*xcol.t());
       for(int j = 0; j < 4; ++j){
           b(j, 0) = b(j, 0) + weight[i]*xx(j, i)/v[i]*(1 - ran[i]*err[i]);
       }
  }

  a = inv(a);
  return(wrap(a*b));
}




// [[Rcpp::export]]
List Boot_Resti(NumericVector weight, NumericVector theta, NumericVector r, String score, int iter = 25, double StSize = 1){
  int dim = theta.size();
  NumericVector a(dim), the(dim);

  the = clone(theta);

  for(int i = 0; i < iter; ++i){
    a = Boot_R_change(weight, the, r, score);
    the = the - a/StSize;
  }

  return List::create(
                      _["change"] = a,
                      _["coef"] = the);
}


// [[Rcpp::export]]
List Boot_Resti_21(NumericVector weight, NumericVector theta, NumericVector r, String score, int iter = 25, double StSize = 1){
  int dim = theta.size();
  NumericVector a(dim), the(dim);

  the = clone(theta);

  for(int i = 0; i < iter; ++i){
    a = Boot_R_change_21(weight, the, r, score);
    the = the - a/StSize;
  }

  return List::create(
                      _["change"] = a,
                      _["coef"] = the);
}




/*QMLE for GARCH(1,1) and GARCH(2,1); compute change at each step*/
NumericVector QMLE_change(NumericVector theta, NumericVector r){
  int n = r.size();
  NumericVector v = vt(theta, r);
  NumericVector y(n), w(n), I1(n), I2(n);
  NumericMatrix x(3, n);

  x = der(theta,r);
  w = 1/pow(v, 2);

  arma::vec rep(n);
  rep = rep.ones();

  y = (pow(r, 2) - v);

  arma::mat a(3, 3), b(3, 1); /*define two matrices*/

  a.zeros();
  b.zeros();

  arma::mat xx = as<arma::mat>(x); /*covert to arma mat*/
  arma::vec xcol(3);

  for(int i = 0; i < n; ++i){
       xcol = xx.col(i);
       a = a + w[i]*(xcol*xcol.t());
       for(int j = 0; j < 3; ++j){
        b(j, 0) = b(j, 0) + w[i]*xx(j, i)*y[i];
       }
  }
       a = inv(a);


  return(wrap(a*b));
}



NumericVector QMLE_change_21(NumericVector theta, NumericVector r){
  int n = r.size();
  NumericVector v = vt_21(theta, r);
  NumericVector y(n), w(n), I1(n), I2(n);
  NumericMatrix x(4, n);

  x = der_21(theta,r);
  w = 1/pow(v, 2);

  arma::vec rep(n);
  rep = rep.ones();

  y = (pow(r, 2) - v);

  arma::mat a(4, 4), b(4, 1); /*define two matrices*/

  a.zeros();
  b.zeros();

  arma::mat xx = as<arma::mat>(x); /*covert to arma mat*/
  arma::vec xcol(4);

  for(int i = 0; i < n; ++i){
       xcol = xx.col(i);
       a = a + w[i]*(xcol*xcol.t());
       for(int j = 0; j < 4; ++j){
        b(j, 0) = b(j, 0) + w[i]*xx(j, i)*y[i];
       }
  }
       a = inv(a);


  return(wrap(a*b));
}


/*QMLE for GARCH(1,1)*/
// [[Rcpp::export]]
List QMLEesti(NumericVector theta, NumericVector r, int iter = 25, double StSize = 1){
  int dim = theta.size();
  NumericVector a(dim), the(dim);

  the = clone(theta);

  for(int i = 0; i < iter; ++i){
    a = QMLE_change(the, r);
    the = the + a/StSize;
  }

  return List::create(
                      _["change"] = a,
                      _["coef"] = the);
}


/*QMLE for GARCH(2,1)*/
// [[Rcpp::export]]
List QMLEesti_21(NumericVector theta, NumericVector r, int iter = 25, double StSize = 1){
  int dim = theta.size();
  NumericVector a(dim), the(dim);

  the = clone(theta);

  for(int i = 0; i < iter; ++i){
    a = QMLE_change_21(the, r);
    the = the + a/StSize;
  }

  return List::create(
                      _["change"] = a,
                      _["coef"] = the);
}




/*LAD for GARCH(1,1)*/


// [[Rcpp::export]]
NumericVector LAD_change(NumericVector theta, NumericVector r){
  int n = r.size();
  NumericVector v = vt(theta, r);
  NumericVector y(n), w(n), I1(n), I2(n);
  NumericMatrix x(3, n);

  x = der(theta,r);
  w = 1/pow(v, 2);

  arma::vec rep(n);
  rep = rep.ones();

  y = pow(v, 0.5)*(abs(r)-pow(v, 0.5));

  arma::mat a(3, 3), b(3, 1); /*define two matrices*/

  a.zeros();
  b.zeros();

  arma::mat xx = as<arma::mat>(x); /*covert to arma mat*/
  arma::vec xcol(3);

  for(int i = 0; i < n; ++i){
       xcol = xx.col(i);
       a = a + w[i]*(xcol*xcol.t());
       for(int j = 0; j < 3; ++j){
        b(j, 0) = b(j, 0) + w[i]*xx(j, i)*y[i];
       }
  }
       a = inv(a);


  return(wrap(a*b));
}


// [[Rcpp::export]]
List LADesti(NumericVector theta, NumericVector r, int iter = 17, double StSize = 1){
  int dim = theta.size();
  NumericVector a(dim), the(dim);

  the = clone(theta);

  for(int i = 0; i < iter; ++i){
    a = LAD_change(the, r);
    the = the + a/StSize;
  }

  return List::create(
                      _["change"] = a,
                      _["coef"] = the);
}








