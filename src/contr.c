/********************************************/
/*        Authors: Hyungwon Choi            */
/*                 Debashis Ghosh           */
/********************************************/
/*       Weighted Contrast Method           */
/*       Ref: J.Wang et.al                  */
/********************************************/
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <Rmath.h>
#include "meta.h"

/***********************/
/*   Main functions    */
/***********************/

void contr(double *exprs, int *labels, int *nd, int *nr, int *nc, int *numperm, double *z, double *p)
{
  /***********************************/
  /*** Variable declaration        ***/
  /***********************************/
  int i;
  ARRAY2 data[*nd];  

  /***********************************/
  /*** initialization of ARRAYS    ***/
  /*** calculates mu and stdev     ***/
  /***********************************/
  init_ARRAYS(exprs,nd,nr,nc,labels,data);  
  
  /************************************/
  /*** perform lowess in each       ***/
  /*** calculate weighted contrasts ***/
  /************************************/
  for(i=0;i<*nd;i++) {
    do_LOWESS(data[i].mean0,data[i].var0,*nr);        
    do_LOWESS(data[i].mean1,data[i].var1,*nr);
  }
  weighted_contrast(data, nd, z, nr); 

  /****************************************/
  /*** p-values referenced to standard  ***/
  /*** normal distribution if permute=0 ***/
  /****************************************/
  if(*numperm == 0) {
    for(i=0;i<*nr;i++) {
      p[i]=(pnorm5(z[i],0.0,1.0,1,0));
      p[i]=(p[i]>0.5 ? 2.0*(1.0-p[i]) : 2.0*p[i]);
    }
  }
  /******************************************/
  /*** p-values referenced to permutation ***/
  /*** based distribution if permute!=0   ***/
  /******************************************/
  else {
    GetRNGstate();   
    permute_pval(data,nd,nr,nc,numperm,z,p);
    PutRNGstate();
  }

  for(i=0;i<*nd;i++) free_array2(&(data[i]));
}



