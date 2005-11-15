/*****************************************/
/*       Authors: Hyungwon Choi          */
/*                Debashis Ghosh         */
/*****************************************/
/*             Data Processing           */
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
/* This is not for POE */
/***********************/

void malloc_array(ARRAY *expr)
{
  static int i, nr, nc;
  nr=expr->nrow;
  nc=expr->ncol;
  /* assert(expr->id=(char **) Calloc(nr, char *)); */
  assert(expr->d=(double **) Calloc(nr, double *));
  assert(expr->label=(int *) Calloc(nc, int));

  /* Initialize */
  memset(expr->label, 0, sizeof(int)*nc);
  for(i=0;i<nc;i++) expr->label[i]=0;
  for(i=0;i<nr;i++) {
    /*  assert(expr->id[i]=(char *) Calloc(MAX_ID, char)); */
    assert(expr->d[i]=(double *) Calloc(nc, double));
  }
}

void malloc_PP(PP *pp, int *nrow, int *ncol) 
{
  static int i, nr, nc;
  nr=*nrow;
  nc=*ncol;
  assert(pp->alpha_t=(double *) Calloc(nc, double));
  assert(pp->mu_g=(double *) Calloc(nr, double));
  assert(pp->kappa_pos_g=(double *) Calloc(nr, double));
  assert(pp->kappa_neg_g=(double *) Calloc(nr, double));
  assert(pp->sigma_g=(double *) Calloc(nr, double));
  assert(pp->pi_pos_g=(double *) Calloc(nr, double));
  assert(pp->pi_neg_g=(double *) Calloc(nr, double));
  assert(pp->poe_mat=(double **) Calloc(nr, double *));
  assert(pp->phat_pos=(double **) Calloc(nr, double *));
  assert(pp->phat_neg=(double **) Calloc(nr, double *));
  for(i=0;i<nr;i++) {
    assert(pp->poe_mat[i]=(double *) Calloc(nc, double));
    assert(pp->phat_pos[i]=(double *) Calloc(nc, double));
    assert(pp->phat_neg[i]=(double *) Calloc(nc, double));
  }
}

void malloc_CH(CH *ch, int *nrow, int *ncol, int *niter)
{
  static int i, j, nr, nc, num;
  nr=*nrow;
  nc=*ncol;
  num=*niter;
  assert(ch->alpha_t=(double **) Calloc(nc,double *));
  assert(ch->mu_g=(double **) Calloc(nr,double *));
  assert(ch->kappa_pos_g=(double **) Calloc(nr,double *));
  assert(ch->kappa_neg_g=(double **) Calloc(nr,double *));
  assert(ch->sigma_g=(double **) Calloc(nr,double *));
  assert(ch->pi_pos_g=(double **) Calloc(nr,double *));
  assert(ch->pi_neg_g=(double **) Calloc(nr,double *));
  assert(ch->poe_mat=(double **) Calloc(nr,double *));
  for(i=0;i<nc;i++) {
    assert(ch->alpha_t[i]=(double *) Calloc(num,double));
  }
  for(i=0;i<nr;i++) {
    assert(ch->mu_g[i]=(double *) Calloc(num,double));
    assert(ch->kappa_pos_g[i]=(double *) Calloc(num,double));
    assert(ch->kappa_neg_g[i]=(double *) Calloc(num,double));
    assert(ch->sigma_g[i]=(double *) Calloc(num,double));
    assert(ch->pi_pos_g[i]=(double *) Calloc(num,double));
    assert(ch->pi_neg_g[i]=(double *) Calloc(num,double));
    assert(ch->poe_mat[i]=(double *) Calloc(nc,double));
  }
  assert(ch->mu=(double *) Calloc(num,double));
  assert(ch->tausqinv=(double *) Calloc(num,double));
  assert(ch->gamma=(double *) Calloc(num,double));
  assert(ch->lambda=(double *) Calloc(num,double));
  assert(ch->pil_pos_mean=(double *) Calloc(num,double));
  assert(ch->pil_neg_mean=(double *) Calloc(num,double));
  assert(ch->pil_pos_prec=(double *) Calloc(num,double));
  assert(ch->pil_neg_prec=(double *) Calloc(num,double));
  assert(ch->kap_pos_rate=(double *) Calloc(num,double));
  assert(ch->kap_neg_rate=(double *) Calloc(num,double));
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      ch->poe_mat[i][j]=0.0;
    }
  }
}

void free_array(ARRAY *expr)
{
  static int i;
  /* free sequentially for all layers */
  for(i=0;i<expr->nrow;i++) {
    Free(expr->d[i]);
  }
  Free(expr->label);
  Free(expr->d);
  /* Free(expr->nrow);
  Free(expr->ncol); */
}

void free_PP(PP *pp, int *nrow)
{
  static int i, nr;
  nr=*nrow;
  Free(pp->alpha_t);
  Free(pp->mu_g);
  Free(pp->kappa_pos_g);
  Free(pp->kappa_neg_g);
  Free(pp->sigma_g);
  Free(pp->pi_pos_g);
  Free(pp->pi_neg_g);
  for(i=0;i<nr;i++) {
    Free(pp->poe_mat[i]);
    Free(pp->phat_pos[i]);
    Free(pp->phat_neg[i]);
  }
  Free(pp->poe_mat);
  Free(pp->phat_pos);
  Free(pp->phat_neg);
}

void free_CH(CH *ch, int *nrow, int *ncol, int *niter)
{
  static int i, num, nr, nc;
  nr=*nrow;
  num=*niter;
  for(i=0;i<nc;i++) Free(ch->alpha_t[i]);
  for(i=0;i<nr;i++) {
    Free(ch->mu_g[i]);
    Free(ch->kappa_pos_g[i]);
    Free(ch->kappa_neg_g[i]);
    Free(ch->sigma_g[i]);
    Free(ch->pi_pos_g[i]);
    Free(ch->pi_neg_g[i]);
    Free(ch->poe_mat[i]);
  }
  Free(ch->alpha_t);
  Free(ch->mu_g);
  Free(ch->kappa_pos_g);
  Free(ch->kappa_neg_g);
  Free(ch->sigma_g);
  Free(ch->pi_pos_g);
  Free(ch->pi_neg_g);
  Free(ch->poe_mat);
  Free(ch->mu);
  Free(ch->tausqinv);
  Free(ch->gamma);
  Free(ch->lambda);
  Free(ch->pil_pos_mean);
  Free(ch->pil_neg_mean);
  Free(ch->pil_pos_prec);
  Free(ch->pil_neg_prec);
  Free(ch->kap_pos_rate);
  Free(ch->kap_neg_rate);
}

/* Still data processing */

void init_ARRAY(double *d, int *nrow, int *ncol, int *label, ARRAY *expr) 
{
  static int i,j;
  expr->nrow=*nrow;
  expr->ncol=*ncol;
  malloc_array(expr);
  for(j=0;j<expr->ncol;j++) expr->label[j]=label[j];
  for(i=0;i<expr->nrow;i++) {
    for(j=0; j<expr->ncol;j++) {
      expr->d[i][j] = d[j*expr->nrow+i];
    }
  }
} 

void init_PP(PP *pp, int *nrow, int *ncol) 
{
  static int i,j,nr,nc;
  nr=*nrow;
  nc=*ncol;
  malloc_PP(pp, nrow, ncol);
  for(i=0;i<nc;i++) pp->alpha_t[i]=0.0;
  for(i=0;i<nr;i++) {
    pp->mu_g[i]=0.0;
    pp->kappa_pos_g[i]=0.0;
    pp->kappa_neg_g[i]=0.0;
    pp->sigma_g[i]=0.0;
    pp->pi_pos_g[i]=0.0;
    pp->pi_neg_g[i]=0.0;
  }
  pp->mu=0.0;
  pp->tausqinv=0.0;
  pp->gamma=0.0;
  pp->lambda=0.0;
  pp->pil_pos_mean=0.0;
  pp->pil_neg_mean=0.0;
  pp->pil_pos_prec=0.0;
  pp->pil_neg_prec=0.0;
  pp->kap_pos_rate=0.0;
  pp->kap_neg_rate=0.0;
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      pp->poe_mat[i][j]=0.0;
      pp->phat_pos[i][j]=0.0;
      pp->phat_neg[i][j]=0.0;
    }
  }
}

void mat2vec(double **mat, double *d, int *nrow, int *ncol) 
{ 
  static int i,j,nr,nc;
  nr=*nrow;
  nc=*ncol;
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      d[j*nr+i]=mat[i][j];
    }
  }
}

void vec2PR(double *vec, PR *pr) 
{
  pr->alpha_mm=vec[0];
  pr->alpha_sd=vec[1];
  pr->mu_mm=vec[2];
  pr->mu_sd=vec[3];
  pr->pi_pos_mm=vec[4];
  pr->pi_pos_sd=vec[5];
  pr->pi_neg_mm=vec[6];
  pr->pi_neg_sd=vec[7];
  pr->kap_pri_rate=vec[8];
  pr->tausqinv_aa=vec[9];
  pr->tausqinv_bb=vec[10];
}

void vec2PP(double *vec, PP *pp, int *nrow, int *ncol)
{
  static int i,j,nr,nc;
  nr=*nrow; 
  nc=*ncol;
  malloc_PP(pp, nrow, ncol);
  for(i=0;i<nc;i++) pp->alpha_t[i]=vec[i];
  for(i=0;i<nr;i++) {
    pp->mu_g[i]=vec[nc+i];
    pp->kappa_pos_g[i]=vec[nc+nr+i];
    pp->kappa_neg_g[i]=vec[nc+(2*nr)+i];
    pp->sigma_g[i]=vec[nc+(3*nr)+i];
    pp->pi_pos_g[i]=vec[nc+(4*nr)+i];
    pp->pi_neg_g[i]=vec[nc+(5*nr)+i];
  }
  pp->mu=vec[nc+(6*nr)];
  pp->tausqinv=vec[nc+(6*nr)+1];
  pp->gamma=vec[nc+(6*nr)+2];
  pp->lambda=vec[nc+(6*nr)+3];
  pp->pil_pos_mean=vec[nc+(6*nr)+4];
  pp->pil_neg_mean=vec[nc+(6*nr)+5];
  pp->pil_pos_prec=vec[nc+(6*nr)+6];
  pp->pil_neg_prec=vec[nc+(6*nr)+7];
  pp->kap_pos_rate=vec[nc+(6*nr)+8];
  pp->kap_neg_rate=vec[nc+(6*nr)+9];
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      pp->poe_mat[i][j]=vec[nc+(6+j)*nr+10+i];
    }
  }
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      pp->phat_pos[i][j]=vec[nc+(6+nc+j)*nr+10+i];
    }
  }
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      pp->phat_neg[i][j]=vec[nc+(6+2*nc+j)*nr+10+i];
    }
  }  
  pp->accept=vec[nc+(6+3*nc)*nr+11];

}

void PP2vec(double *vec, PP *res, int *nrow, int *ncol)
{
  static int i,j,nr,nc;
  nr=*nrow;
  nc=*ncol;
  for(i=0;i<nc;i++) vec[i]=res->alpha_t[i];
  for(i=0;i<nr;i++) {
    vec[nc+i]=res->mu_g[i];
    vec[nc+nr+i]=res->kappa_pos_g[i];
    vec[nc+(2*nr)+i]=res->kappa_neg_g[i];
    vec[nc+(3*nr)+i]=res->sigma_g[i];
    vec[nc+(4*nr)+i]=res->pi_pos_g[i];
    vec[nc+(5*nr)+i]=res->pi_neg_g[i];
  }
  vec[nc+(6*nr)]=res->mu;
  vec[nc+(6*nr)+1]=res->tausqinv;
  vec[nc+(6*nr)+2]=res->gamma;
  vec[nc+(6*nr)+3]=res->lambda;
  vec[nc+(6*nr)+4]=res->pil_pos_mean;
  vec[nc+(6*nr)+5]=res->pil_neg_mean;
  vec[nc+(6*nr)+6]=res->pil_pos_prec;
  vec[nc+(6*nr)+7]=res->pil_neg_prec;
  vec[nc+(6*nr)+8]=res->kap_pos_rate;
  vec[nc+(6*nr)+9]=res->kap_neg_rate;
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      vec[nc+(6*nr)+10+(j*nr)+i]=res->poe_mat[i][j];
    }
  }
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      vec[nc+(nc+j+6)*nr+10+i]=res->phat_pos[i][j];
    }
  }
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      vec[nc+(2*nc+j+6)*nr+10+i]=res->phat_neg[i][j];
    } 
  }
  vec[nc+(3*nc+6)*nr+11]=res->accept;

}

/**************************************/
/* summarizing results from CH struct */
/**************************************/

void update_CH(CH *ch, PP *pp, int iter, int *numiter, int *nrow, int *ncol)
{
  static int i,j,n,nr,nc;
  n=iter;
  nr=*nrow;
  nc=*ncol;
  for(j=0;j<nc;j++) {
    ch->alpha_t[j][n]=pp->alpha_t[j];
  }
  for(i=0;i<nr;i++) {
    ch->mu_g[i][n]=pp->mu_g[i];
    ch->kappa_pos_g[i][n]=pp->kappa_pos_g[i];
    ch->kappa_neg_g[i][n]=pp->kappa_neg_g[i];
    ch->sigma_g[i][n]=pp->sigma_g[i];
    ch->pi_pos_g[i][n]=pp->pi_pos_g[i];
    ch->pi_neg_g[i][n]=pp->pi_neg_g[i];
  }
  ch->mu[n]=pp->mu;
  ch->tausqinv[n]=pp->tausqinv;
  ch->gamma[n]=pp->gamma;
  ch->lambda[n]=pp->lambda;
  ch->pil_pos_mean[n]=pp->pil_pos_mean;
  ch->pil_neg_mean[n]=pp->pil_neg_mean;
  ch->pil_pos_prec[n]=pp->pil_pos_prec;
  ch->pil_neg_mean[n]=pp->pil_neg_prec;
  ch->kap_pos_rate[n]=pp->kap_pos_rate;
  ch->kap_neg_rate[n]=pp->kap_neg_rate;
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      ch->poe_mat[i][j]+=pp->poe_mat[i][j]/((double) *numiter);
    }
  }
  ch->accept+=pp->accept/((double) *numiter);
}
 
void median_CH(CH *ch, PP *res, int len, int *nrow, int *ncol)
{
  static int nr, nc, i, j;
  nr=*nrow;
  nc=*ncol;
  for(j=0;j<nc;j++) {
    res->alpha_t[j]=get_median(ch->alpha_t[j],len);
  }
  for(i=0;i<nr;i++) {
    res->mu_g[i]=get_median(ch->mu_g[i],len);
    res->kappa_pos_g[i]=get_median(ch->kappa_pos_g[i],len);
    res->kappa_neg_g[i]=get_median(ch->kappa_neg_g[i],len);
    res->sigma_g[i]=get_median(ch->sigma_g[i],len);
    res->pi_pos_g[i]=get_median(ch->pi_pos_g[i],len);
    res->pi_neg_g[i]=get_median(ch->pi_neg_g[i],len);
  }
  res->mu=get_median(ch->mu,len);
  res->tausqinv=get_median(ch->tausqinv,len);
  res->gamma=get_median(ch->gamma,len);
  res->lambda=get_median(ch->lambda,len);
  res->pil_pos_mean=get_median(ch->pil_pos_mean,len);
  res->pil_neg_mean=get_median(ch->pil_neg_mean,len);
  res->pil_pos_prec=get_median(ch->pil_pos_prec,len);
  res->pil_neg_prec=get_median(ch->pil_neg_prec,len);
  res->kap_pos_rate=get_median(ch->kap_pos_rate,len);
  res->kap_neg_rate=get_median(ch->kap_neg_rate,len); 
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      res->poe_mat[i][j]=ch->poe_mat[i][j];
    }
  }
  res->accept=ch->accept;
 
}

double get_median(double *base, int len)
{
  double *new_vec;
  double med;
  int i, pk;
  assert(new_vec=(double *) Calloc(len,double));  
  for(i=0;i<len;i++) {
    new_vec[i]=base[i]; 
  }
  /* pk=(len%2==0 ? len/2 : (len+1)/2); */
  if(len==1) {
    med=base[0];
    Free(new_vec);
    return med;
  }
  else if(len%2==0) {
    R_rsort(new_vec,len);
    pk=(len-2)/2;
    med=(new_vec[pk]+new_vec[pk+1])/2.0;
    Free(new_vec);
    return med;
  }
  else {
    R_rsort(new_vec,len);
    pk=(len-1)/2;
    med=new_vec[pk];
    Free(new_vec);
    return med;
  }
}



/**************************************************/
/* Functions needed in Weighted Contrasts methods */
/**************************************************/
/* See dproc2.c                                   */
/**************************************************/
