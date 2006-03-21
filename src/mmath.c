/*****************************************/
/*       Authors: Hyungwon Choi          */
/*                Debashis Ghosh         */
/*****************************************/
/*             Math Functions            */
/*****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>
#include <Rmath.h>
#include "meta.h"
#define logit(p) log(p)-log(1-p)

double vec_sum(const double *vec, int len) {
  static int i;
  static double res;
  res=vec[0];
  for(i=1;i<len;i++) res+=vec[i];
  return res;
}

double vec_max(const double *vec, int len) {
  static int i;
  static double res;
  res=vec[0];
  for(i=1;i<len;i++) {
    if(res<vec[i]) res=vec[i];
  }
  return res;
}

double vec_min(const double *vec, int len) {
  static int i;
  static double res;
  res=vec[0];
  for(i=1;i<len;i++) {
    if(res>vec[i]) res=vec[i];
  }
  return res;
}

double vec_mean(const double *vec, int len) {
  static double tmp=0.0;
  static int i;
  for(i=0;i<len;i++) tmp+=vec[i];
  tmp=tmp/((double) len);
  return tmp;
}

double vec_var(const double *vec, int len) {
  static double mean=0.0;
  static double var=0.0;
  static int i;
  for(i=0;i<len;i++) mean+=vec[i];
  mean=mean/((double) len);
  for(i=0;i<len;i++) var+=pow((vec[i]-mean),2);
  var/=((double) (len-1));
  var=sqrt(var);
  return var;
}

double xlogy(const double x, const double y) 
{
  double y2;
  y2=(((y<=0.0) & (y > -0.001)) ? 1e-8 : y);
  return x*log(y2);
}

double log_posterior_gamma(const double a, const double b, const double *sigma_g, const int len)
{
  static double lp=0.0;
  static int i;
  if(a < 0.0) return _NA_VAL_;
  else {
    for(i=0;i<len;i++) {
      lp+=-lgammafn(a)+a*log(b)+(a-1)*log(1.0/pow(sigma_g[i],2.0));
    }
    return lp;
  }
}

double log_posterior_kappa(const double kappa, const double *xxx, const int len, const double sigma_g, const double pi_g, const double kr)
{
  static int i;
  double lp, lik, sum=0.0;
  lik=1.0;
  if(len==0) {
    lp=log(kr)-kr*kappa;
    return lp;
  }
  if(kappa<0) return _NA_VAL_;
  for(i=0;i<len;i++) sum+=((int) (xxx[i]<=kappa));
  if(sum==len) {
    for(i=0;i<len;i++) {
      lik*=pi_g*.5*(1.0/kappa)+(1-pi_g)*dnorm4(xxx[i],0,sigma_g,0);
    }
  }
  else {
    sum=0.0;
    for(i=0;i<len;i++) sum+=((int) (xxx[i]>kappa));
    if(sum==len) {
      for(i=0;i<len;i++) {
        lik*=(1-pi_g)*dnorm4(xxx[i],0,sigma_g,0);
      }
    }
    else {
      for(i=0;i<len;i++) {
        if(xxx[i]<=kappa) {
          lik*=pi_g*.5*(1.0/kappa)+(1-pi_g)*dnorm4(xxx[i],0,sigma_g,0);   
	}
        else{
          lik*=(1-pi_g)*dnorm4(xxx[i],0,sigma_g,0);
        }
      }
    }
  }
  lp=log(lik)+log(kr)-kr*kappa;
  return lp;
}

void calcor(const double *x, const double *y, const int len, double *num)
{
  static int i;
  static double correl, mx, my;
  static double qx, qy, xy;
  static double t0, t1, t2;
  correl=0.0;
  t0=0.0;
  t1=0.0;
  t2=0.0;
  mx=0.0;
  my=0.0;
  qx=0.0;
  qy=0.0;
  xy=0.0;
  for(i=0;i<len;i++) {
    mx+=x[i];
    my+=y[i];
    qx+=x[i]*x[i];
    qy+=y[i]*y[i];
    xy+=x[i]*y[i];
  }
  t0 = ((double) xy) - (mx*my/((double)len));
  t1 = ((double) qx) - ((mx*mx)/((double)len));
  t2 = ((double) qy) - ((my*my)/((double)len));
  correl = t0 / (sqrt(t1) * sqrt(t2));
  *num = correl;
}


