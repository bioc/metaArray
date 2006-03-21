/*****************************************/
/*       Authors: Hyungwon Choi          */
/*                Debashis Ghosh         */
/*****************************************/
/* Data Processing for Weighted Contrast */
/*****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>
#include <Rmath.h>
#include "meta.h"

/***********************/
/* Memory allocation   */
/***********************/

void malloc_array2(ARRAY2 *expr)
{
  static int i, nr, nc;
  nr=expr->nrow;
  nc=expr->ncol;
  assert(expr->d=(double **) Calloc(nr, double *));
  assert(expr->label=(int *) Calloc(nc, int));
  assert(expr->mean0=(double *) Calloc(nr, double));
  assert(expr->var0=(double *) Calloc(nr, double));
  assert(expr->mean1=(double *) Calloc(nr, double));
  assert(expr->var1=(double *) Calloc(nr, double));
  assert(expr->mean_diff=(double *) Calloc(nr, double));
  assert(expr->var_sum=(double *) Calloc(nr, double));
  /* Initialize */
  /* memset(expr->label, 0, sizeof(int)*nc); */
  for(i=0;i<nc;i++) expr->label[i]=0;
  for(i=0;i<nr;i++) {
    assert(expr->d[i]=(double *) Calloc(nc, double));
  }
}


void free_array2(ARRAY2 *expr)
{
  static int i;
  /* free sequentially for all layers */
  for(i=0;i<expr->nrow;i++) {
    Free(expr->d[i]);
  }
  Free(expr->label);
  Free(expr->mean0);
  Free(expr->var0);
  Free(expr->mean1);
  Free(expr->var1);
  Free(expr->mean_diff);
  Free(expr->var_sum);
  Free(expr->d);
  /* Free(expr->nrow);
  Free(expr->ncol); */
}

void get_meanvar(ARRAY2 *expr)
{
  int i,j,n0,n1;
  n0=0; n1=0;
  for(j=0;j<expr->ncol;j++) {
    if(expr->label[j]==0) n0++;
    if(expr->label[j]==1) n1++;
  }
  for(i=0;i<expr->nrow;i++) {
    expr->mean0[i]=0.0;
    expr->var0[i]=0.0;
    expr->mean1[i]=0.0;
    expr->var1[i]=0.0;
  }
  for(i=0;i<expr->nrow;i++) {
    for(j=0;j<expr->ncol;j++) {
      if(expr->label[j]==0) expr->mean0[i] += expr->d[i][j];
      if(expr->label[j]==1) expr->mean1[i] += expr->d[i][j];
    }
    expr->mean0[i] /= ((double) n0);
    expr->mean1[i] /= ((double) n1);
    /* expr->mean_diff[i] = expr->mean1[i] - expr->mean0[i]; */
    for(j=0;j<expr->ncol;j++) {
      if(expr->label[j]==0) expr->var0[i] += pow((expr->d[i][j]-expr->mean0[i]),2);
      if(expr->label[j]==1) expr->var1[i] += pow((expr->d[i][j]-expr->mean1[i]),2);
    }
    expr->var0[i] /= ((double) (n0 - 1));
    expr->var1[i] /= ((double) (n1 - 1));
    /* expr->var_sum[i] = expr->var0[i]/((double) n0) + expr->var1[i]/((double) n1);  */
  }
}

void init_ARRAY2(double *d, int *nrow, int *ncol, int *label, ARRAY2 *expr) 
{
  static int i,j;
  expr->nrow=*nrow;
  expr->ncol=*ncol;
  malloc_array2(expr);
  for(j=0;j<expr->ncol;j++) {
    expr->label[j]=label[j];
  }
  for(i=0;i<expr->nrow;i++) {
    for(j=0;j<expr->ncol;j++) {
      expr->d[i][j] = d[j*expr->nrow+i];
    }
  }
  get_meanvar(expr);
}

void init_ARRAYS(double *exprs, int *ndata, int *nrow, int *ncol, int *labels, ARRAY2 data[])
{
  static int i,j,k,cum1,cum2;
  static int *cl;
  static double *expr;
  cum1=0;
  cum2=0;
  for(i=0;i<*ndata;i++) {
    expr = (double *) Calloc((*nrow)*(ncol[i]),double);
    cl = (int *) Calloc(ncol[i],int);
    for(j=0;j<ncol[i];j++) {
      for(k=0;k<*nrow;k++) expr[j*(*nrow)+k] = exprs[cum1+j*(*nrow)+k];
      cl[j] = labels[cum2+j];
    }
    init_ARRAY2(expr,nrow,&(ncol[i]),cl,&(data[i]));
    Free(expr);
    Free(cl);
    cum2+=ncol[i];
    cum1=cum2*(*nrow);
  }
}

void do_LOWESS(double *x, double *y, int len) 
{
  static double *tx, *ty /* *tmp  */ ;
  static double *ys, *rw, *res;
  static int *ind;
  static int i,j, nsteps,k;
  static double delta, f /* max, min  */;
  f=2.0/3.0;
  delta = (vec_max(x,len)-vec_min(x,len))*0.01;
  nsteps=3;
  assert(ind = (int *) Calloc(len,int));
  assert(tx = (double *) Calloc(len,double));
  assert(ty = (double *) Calloc(len,double));
  assert(ys = (double *) Calloc(len,double));
  assert(rw = (double *) Calloc(len,double));
  assert(res = (double *) Calloc(len,double));
  for(i=0;i<len;i++) {
    ind[i]=i;
    tx[i]=x[i];  
    ty[i]=y[i];
  }
  memset(ys, 0.0, sizeof(double)*len);
  memset(rw, 0.0, sizeof(double)*len);
  memset(res, 0.0, sizeof(double)*len); 
  rsort_with_index(tx,ind,len);    
  for(i=0;i<len;i++) ty[i]=y[ind[i]];    
  lowess(tx,ty,&len,&f,&nsteps,&delta,ys,rw,res);
  for(i=0;i<len;i++) {
    k=0;
    for(j=0;(i<len && k==0);j++) {
      if(x[i]==tx[j]) {
        y[i]=ys[j];
        k=1;
      }
    }
  }
  Free(ind);
  Free(tx);
  Free(ty);
  Free(ys);
  Free(rw);
  Free(res);
}

void weighted_contrast(ARRAY2 data[], int *nd, double *z, int *nrow) {
  static int i,j,k;
  static double *denom;
  static double diff, va;
  int n0,n1;
  assert(denom=(double *) Calloc(*nrow,double));
  for(i=0;i<*nrow;i++) {
    z[i] = 0.0;
    denom[i] = 0.0;
  }
  for(i=0;i<*nrow;i++) {
    for(j=0;j<*nd;j++) {
      n0=0; n1=0;
      for(k=0;k<data[j].ncol;k++) {
        if(data[j].label[k]==0) n0++;
        if(data[j].label[k]==1) n1++;
      }
      diff = data[j].mean1[i] - data[j].mean0[i];
      va = data[j].var0[i]/((double) n0) + data[j].var1[i]/((double) n1);
      z[i] += diff / va;
      denom[i] += ((double) 1.0) / (va);     
    }
    z[i] /= sqrt(denom[i]);
  }
  Free(denom);
}


void permute_pval(ARRAY2 data[], int *nd, int *nr, int *nc, int *numperm, double *z,double *p) {
  double **permu;
  double *zz;
  int **cl;
  int i,j,k;
  int tmp;

  assert(zz = (double *) Calloc(*nr,double)); 
  assert(permu = (double **) Calloc(*numperm,double *));
  for(i=0;i<*numperm;i++) {
    assert(permu[i] = (double *) Calloc(*nr,double));
  }
  assert(cl = (int **) Calloc(*nd,int *));
  for(i=0;i<*nd;i++) {    
    assert(cl[i] = (int *) Calloc(nc[i],int));
  }
  for(i=0;i<*nd;i++) {
    for(j=0;j<nc[i];j++) {
      cl[i][j] = data[i].label[j];
    }
  }
  
  /**************************************************************/
  /*** Feed randomly permuted class labels in each meta array ***/
  /*** and calculate contrast scores again and thus construct ***/
  /*** nonparametric distribution of weighted scores per gene ***/
  /*** Finally, count the number of permuted samples greater  ***/
  /*** than observed value in absolute value.                 ***/
  /**************************************************************/

  for(k=0;k<*numperm;k++) {
    for(i=0;i<*nd;i++) {
      perm(cl[i],data[i].label,nc[i]);
      get_meanvar(&(data[i]));
    }     
    for(i=0;i<*nr;i++) zz[i]=0.0;
    weighted_contrast(data, nd, zz, nr); 
    for(i=0;i<*nr;i++) {
      permu[k][i] = zz[i];
    }
  }
 
  for(i=0;i<*nr;i++) {
    for(j=0;j<*numperm;j++) {
      tmp = permu[j][i] > z[i] ? 1 : 0 ;
      p[i] += ((double) tmp) /((double) *numperm) ;
    }
  } 
  for(i=0;i<*nr;i++) p[i]=(p[i]>0.5 ? 2.0*(1.0-p[i]) : 2.0*p[i]);

  for(i=0;i<*numperm;i++) Free(permu[i]);
  Free(permu);
  Free(zz);

}



