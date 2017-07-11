#include <Rcpp.h>
#include <cmath>  
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector p01CPP(NumericVector x, bool l) {
  int n = x.size();
  NumericVector z(n); 
  NumericVector prob(n);
  for(int i = 0; i < n; ++i) {
    z[i] = cos(x[i]);
    if (z[i]<1){
      if (fabs(z[i])>pow(10.0,-8.0)){
        prob[i] = -1/(2*M_PI*pow(z[i],2.0))*(z[i]+log(1-z[i]));
      } else {
        prob[i] = 1/(4*M_PI);}
    } else {
      prob[i] = pow(10.0,-16.0);
    }
  }
  if (l) {
    prob = log(prob);
  }
  return(prob);
}

// [[Rcpp::export]]
NumericVector p010CPP(NumericVector x,NumericVector y, bool l) {
  NumericVector z = p01CPP(y,l);
  int n = x.size();
  int m = y.size();
  NumericMatrix phi(n, m);
  NumericMatrix joint(n, m);
  NumericMatrix prob(n,m);
  bool bool_valid[n][m]; 
  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      phi(i,j) = x[i]+y[j];
      bool_valid[i][j] = (fabs(phi(i,j))<M_PI)&&(x[i]*y[j]>0);
      if (fabs(phi(i,j))>pow(10.0,-16.0)){
        joint(i,j) = 1/(2*M_PI)*fabs(sin(phi(i,j)))*bool_valid[i][j]/pow(sin(phi(i,j))*cos(x[i])+sin(x[i])*(1-cos(phi(i,j))),2.0); 
      } else {
        joint(i,j) = 1/(2*M_PI)*bool_valid[i][j];
      }
      if (l) {
        prob(i,j) = log(joint(i,j)) - z[j];
      } else{
        prob(i,j) = joint(i,j)/z[j];
      }
    }
  }
  return prob;
}

double my_convolution(double t, double x, double y){
  NumericVector temp_vec1 = p01CPP(wrap(t),0);
  NumericVector temp_vec2 = p010CPP(wrap(x-t),wrap(y),0);
  double out = temp_vec1[0]*temp_vec2[0]; 
  return out;
}

double integral(double x, double y, double a, double b, int n) {
  // this specifically integrates the function my_convolution
  // rather than passing a general function to integrate
  double step = (b - a) / n;  // width of each small rectangle
  double area = 0.0;  // signed area
  for (int i = 0; i < n; i ++) {
    double xi = a + (i + 0.5) * step;
    area += my_convolution(xi, x, y) * step; // sum up each small rectangle
  }
  return area;
}

// [[Rcpp::export]]
NumericVector p011CPP(NumericVector x,NumericVector y, bool l) {
  int n = x.size();
  int m = y.size();
  NumericMatrix prob(n,m);  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      prob(i,j) = integral(x[i], y[j], x[i]-M_PI, x[i]+M_PI, 100);
      if (l) {
        prob(i,j) = log(prob(i,j));
      }
    }
  }
  return wrap(prob); 
}

// noisy emission probabilities calculated via convolution

double my_convolution_q01(double t, double x, double sigma){
  double temp_vec1 = 1/sqrt(2*M_PI*sigma)*exp(-pow(t,2.0)/(2*sigma));
  NumericVector temp_vec2 = p01CPP(wrap(x-t),0);
  double out = temp_vec1*temp_vec2[0]; 
  return out;
}

double integral_q01(double x, double sigma, double a, double b, int n) {
  // this specifically integrates the function my_convolution
  // rather than passing a general function to integrate
  double step = (b - a) / n;  // width of each small rectangle
  double area = 0.0;  // signed area
  for (int i = 0; i < n; i ++) {
    double xi = a + (i + 0.5) * step;
    area += my_convolution_q01(xi, x, sigma) * step; // sum up each small rectangle
  }
  return area;
}

// [[Rcpp::export]]
NumericVector q01CPP(NumericVector x, double sigma, bool l) {
  int n = x.size();
  NumericVector prob(n);  
  for(int i = 0; i < n; ++i) {
      prob[i] = integral_q01(x[i], sigma, x[i]-M_PI, x[i]+M_PI, 100);
      if (l) {
        prob[i] = log(prob[i]);
      }
  }
  return prob; 
}

// noisy emission probabilities calculated via convolution

double my_convolution_q010(double t, double x, double y, double sigma){
  double temp_vec1 = 1/sqrt(2*M_PI*sigma)*exp(-pow(t,2.0)/(2*sigma));
  NumericVector temp_vec2 = p010CPP(wrap(x-t),wrap(y),0);
  double out = temp_vec1*temp_vec2[0]; 
  return out;
}

double integral_q010(double x, double y, double sigma, double a, double b, int n) {
  // this specifically integrates the function my_convolution
  // rather than passing a general function to integrate
  double step = (b - a) / n;  // width of each small rectangle
  double area = 0.0;  // signed area
  for (int i = 0; i < n; i ++) {
    double xi = a + (i + 0.5) * step;
    area += my_convolution_q010(xi, x, y, sigma) * step; // sum up each small rectangle
  }
  return area;
}

// [[Rcpp::export]]
NumericVector q010CPP(NumericVector x, NumericVector y, double sigma, bool l) {
  int n = x.size();
  int m = y.size();
  NumericMatrix prob(n,m);  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      prob(i,j) = integral_q010(x[i], y[j], sigma, x[i]-M_PI, x[i]+M_PI, 100);
      if (l) {
        prob(i,j) = log(prob(i,j));
      }
    }
  }
  return wrap(prob); 
}

double my_convolution_q011(double t, double x, double y, double sigma){
  double temp_vec1 = 1/sqrt(2*M_PI*sigma)*exp(-pow(t,2.0)/(2*sigma));
  NumericVector temp_vec2 = p011CPP(wrap(x-t),wrap(y),0);
  double out = temp_vec1*temp_vec2[0]; 
  return out;
}

double integral_q011(double x, double y, double sigma, double a, double b, int n) {
  // this specifically integrates the function my_convolution
  // rather than passing a general function to integrate
  double step = (b - a) / n;  // width of each small rectangle
  double area = 0.0;  // signed area
  for (int i = 0; i < n; i ++) {
    double xi = a + (i + 0.5) * step;
    area += my_convolution_q011(xi, x, y, sigma) * step; // sum up each small rectangle
  }
  return area;
}

// [[Rcpp::export]]
NumericVector q011CPP(NumericVector x, NumericVector y, double sigma, bool l) {
  int n = x.size();
  int m = y.size();
  int nk = 4*M_PI/sqrt(sigma); //number of integration pts
  NumericMatrix prob(n,m);  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      prob(i,j) = integral_q011(x[i], y[j], sigma, x[i]-M_PI, x[i]+M_PI, nk);
      if (l) {
        prob(i,j) = log(prob(i,j));
      }
    }
  }
  return wrap(prob); 
}

