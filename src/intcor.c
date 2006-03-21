/********************************************/
/*        Authors: Hyungwon Choi            */
/*                 Debashis Ghosh           */
/********************************************/
/*     Integrative Correlation Method       */
/*       Ref: Parmigiani et.al              */
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

void intcor(double *exprs, int *labels, int *nd, int *nr, int *nc, double *correl, double *paircor)
{
  /***********************************/
  /*** Variable declaration        ***/
  /***********************************/
  int i,j,k,ii,jj,ct,id,comb;
  ARRAY2 data[*nd];  
  double **vec; 
  double *tmp1, *tmp2;
  double tmp;
  ct=0;

  /***********************************/
  /*** initialization of ARRAYS    ***/
  /*** calculates mu and stdev     ***/
  /***********************************/
  init_ARRAYS(exprs,nd,nr,nc,labels,data);  
  assert(vec = (double **) Calloc(*nd, double *));
  for(i=0;i<*nd;i++) {
    assert(vec[i] = (double *) Calloc(*nr-1, double));
  }  
  Rprintf("%s", "Gene-specific Integrative Correlations\n");

  /***********************************/
  /*** Begin calculation           ***/
  /***********************************/
  for(i=0;i<*nr;i++) {
    for(ii=0;ii<*nd;ii++) {
      for(jj=0;jj<*nr-1;jj++) {
        vec[ii][jj] = 0.0;
      }
    }
    for(j=0;j<*nd;j++) {
      assert(tmp1=(double *) Calloc(data[j].ncol, double)); 
      assert(tmp2=(double *) Calloc(data[j].ncol, double));
      ct=0; 
      id=0;
      for(k=0;k<data[j].ncol;k++) tmp1[k]=data[j].d[i][k];
      while(ct < *nr) {
        if(ct == i) id=1;
        else {
          for(k=0;k<data[j].ncol;k++) {
            tmp2[k]=data[j].d[ct][k];
          }
          calcor(tmp1, tmp2, data[j].ncol,&(vec[j][ct-id]));
        }
        ct++;
      }
      Free(tmp1);
      Free(tmp2);
    }
    comb = (*nd) * (*nd-1) / 2;    
    tmp=0.0; 
    correl[i]=0.0;
    ct = 0;
    for(ii=0;ii<(*nd)-1;ii++) {
      for(jj=ii+1;jj<*nd ; jj++) {
        calcor(vec[ii],vec[jj],*nr-1, &tmp);
        paircor[ct*(*nr)+i] = tmp;
        tmp = tmp/((double) comb);
        correl[i]+=tmp;
        ct++;
      }
    if((i%100==0) & (i>0)) Rprintf("%i%s", i, " ");
    if((i%1000==0) & (i>0)) Rprintf("%s", "\n");
    }
  }
  for(i=0;i<*nd;i++) Free(vec[i]);
  Free(vec); 
  for(i=0;i<*nd;i++) free_array2(&(data[i]));
}


