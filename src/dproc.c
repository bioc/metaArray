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
  int i, nr, nc;
  nr=expr->nrow;
  nc=expr->ncol;
  /* expr->id=(char **) Calloc(nr, char *); */
  expr->d=(double **) Calloc(nr, double *);
  expr->label=(int *) Calloc(nc, int);

  /* Initialize */
  memset(expr->label, 0, sizeof(int)*nc);
  for(i=0;i<nc;i++) expr->label[i]=0;
  for(i=0;i<nr;i++) {
    /*  expr->id[i]=(char *) Calloc(MAX_ID, char); */
    expr->d[i]=(double *) Calloc(nc, double);
  }
}

void malloc_PP(PP *pp, int *nrow, int *ncol) 
{
  int i, nr, nc;
  nr=*nrow;
  nc=*ncol;
  pp->alpha_t=(double *) Calloc(nc, double);
  pp->mu_g=(double *) Calloc(nr, double);
  pp->kappa_pos_g=(double *) Calloc(nr, double);
  pp->kappa_neg_g=(double *) Calloc(nr, double);
  pp->sigma_g=(double *) Calloc(nr, double);
  pp->pi_pos_g=(double *) Calloc(nr, double);
  pp->pi_neg_g=(double *) Calloc(nr, double);
  pp->poe_mat=(double **) Calloc(nr, double *);
  pp->phat_pos=(double **) Calloc(nr, double *);
  pp->phat_neg=(double **) Calloc(nr, double *);
  for(i=0;i<nr;i++) {
    pp->poe_mat[i]=(double *) Calloc(nc, double);
    pp->phat_pos[i]=(double *) Calloc(nc, double);
    pp->phat_neg[i]=(double *) Calloc(nc, double);
  }
}

/*
void malloc_CH(CH *ch, int *nrow, int *ncol, int *niter)
{
  int i, j, nr, nc, num;
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
  ch->accept = 0.0;
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      ch->poe_mat[i][j]=0.0;
    }
  }
}
*/

void malloc_CH(CH *ch, int *nrow, int *ncol, int *niter)
{
  int i, j, nr, nc, num;
  nr=*nrow;
  nc=*ncol;
  num=*niter;
  ch->alpha_t=(double *) Calloc(nc,double);
  ch->mu_g=(double *) Calloc(nr,double );
  ch->kappa_pos_g=(double *) Calloc(nr,double );
  ch->kappa_neg_g=(double *) Calloc(nr,double );
  ch->sigma_g=(double *) Calloc(nr,double );
  ch->pi_pos_g=(double *) Calloc(nr,double );
  ch->pi_neg_g=(double *) Calloc(nr,double );
  ch->poe_mat=(double **) Calloc(nr,double *);
  for(i=0;i<nr;i++) {
    ch->poe_mat[i]=(double *) Calloc(nc,double);
  }

  ch->accept = 0.0;
  for(j=0;j<nc;j++) ch->alpha_t[j] = 0.0;
  for(i=0;i<nr;i++) {
    ch->mu_g[i] = 0.0;
    ch->kappa_pos_g[i] = 0.0;
    ch->kappa_neg_g[i] = 0.0;
    ch->sigma_g[i] = 0.0;
    ch->pi_pos_g[i] = 0.0;
    ch->pi_neg_g[i] = 0.0;
  }
  ch->mu = 0.0;
  ch->tausqinv = 0.0;
  ch->gamma = 0.0;
  ch->lambda = 0.0;
  ch->pil_pos_mean = 0.0;
  ch->pil_neg_mean = 0.0;
  ch->pil_pos_prec = 0.0;
  ch->pil_neg_prec = 0.0;
  ch->kap_pos_rate = 0.0;
  ch->kap_neg_rate = 0.0;

  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      ch->poe_mat[i][j]=0.0;
    }
  }
}


void free_array(ARRAY *expr)
{
  int i;
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
  int i, nr;
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
  int num, nr, nc;
  nr=*nrow;
  num=*niter;
  Free(ch->alpha_t);
    Free(ch->mu_g);
    Free(ch->kappa_pos_g);
    Free(ch->kappa_neg_g);
    Free(ch->sigma_g);
    Free(ch->pi_pos_g);
    Free(ch->pi_neg_g);
    Free(ch->poe_mat);

}

/* Still data processing */

void init_ARRAY(double *d, int *nrow, int *ncol, int *label, ARRAY *expr) 
{
  int i,j;
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
  int i,j,nr,nc;
  nr=*nrow;
  nc=*ncol;
  malloc_PP(pp, nrow, ncol);
  for(i=0;i<nc;i++) pp->alpha_t[i]=0.0;
  for(i=0;i<nr;i++) {
    pp->mu_g[i]=0.0;
    pp->kappa_pos_g[i]=2.0;
    pp->kappa_neg_g[i]=2.0;
    pp->sigma_g[i]=0.0;
    pp->pi_pos_g[i]=0.2;
    pp->pi_neg_g[i]=0.2;
  }
  pp->mu=0.0;
  pp->tausqinv=1.0;
  pp->gamma=1.0;
  pp->lambda=1.0;
  pp->pil_pos_mean=0.0;
  pp->pil_neg_mean=0.0;
  pp->pil_pos_prec=0.0;
  pp->pil_neg_prec=0.0;
  pp->kap_pos_rate=1.0;
  pp->kap_neg_rate=1.0;
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      pp->poe_mat[i][j]=0.0;
      pp->phat_pos[i][j]=0.2;
      pp->phat_neg[i][j]=0.2;
    }
  }
}

void mat2vec(double **mat, double *d, int *nrow, int *ncol) 
{ 
  int i,j,nr,nc;
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
  int i,j,nr,nc;
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
  int i,j,nr,nc;
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
  int i,j,n,nr,nc;
  n = (*numiter) / 10;
  nr=*nrow;
  nc=*ncol;
  for(j=0;j<nc;j++) {
    ch->alpha_t[j]+=pp->alpha_t[j] / ((double) n);
  }
  for(i=0;i<nr;i++) {
    ch->mu_g[i]+=pp->mu_g[i] / ((double) n);
    ch->kappa_pos_g[i]+=pp->kappa_pos_g[i] / ((double) n);
    ch->kappa_neg_g[i]+=pp->kappa_neg_g[i] / ((double) n);
    ch->sigma_g[i]+=pp->sigma_g[i] / ((double) n);
    ch->pi_pos_g[i]+=pp->pi_pos_g[i] / ((double) n);
    ch->pi_neg_g[i]+=pp->pi_neg_g[i] / ((double) n);
  }
  ch->mu += pp->mu / ((double) n);
  ch->tausqinv += pp->tausqinv / ((double) n);
  ch->gamma += pp->gamma / ((double) n);
  ch->lambda += pp->lambda / ((double) n);
  ch->pil_pos_mean += pp->pil_pos_mean / ((double) n);
  ch->pil_neg_mean += pp->pil_neg_mean / ((double) n);
  ch->pil_pos_prec += pp->pil_pos_prec / ((double) n);
  ch->pil_neg_mean += pp->pil_neg_prec / ((double) n);
  ch->kap_pos_rate += pp->kap_pos_rate / ((double) n);
  ch->kap_neg_rate += pp->kap_neg_rate / ((double) n);
  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      ch->poe_mat[i][j] += pp->poe_mat[i][j] / ((double) n);
    }
  }
  ch->accept += pp->accept / ((double) n);
}
 
void median_CH(CH *ch, PP *res, int len, int *nrow, int *ncol)
{

  /* removed get_median */
  int nr, nc, i, j;
  nr=*nrow;
  nc=*ncol;
  for(j=0;j<nc;j++) {
    res->alpha_t[j]=(ch->alpha_t[j]);
  }
  for(i=0;i<nr;i++) {
    res->mu_g[i]=(ch->mu_g[i]);
    res->kappa_pos_g[i]=(ch->kappa_pos_g[i]);
    res->kappa_neg_g[i]=(ch->kappa_neg_g[i]);
    res->sigma_g[i]=(ch->sigma_g[i]);
    res->pi_pos_g[i]=(ch->pi_pos_g[i]);
    res->pi_neg_g[i]=(ch->pi_neg_g[i]);
  }
  res->mu=(ch->mu);
  res->tausqinv=(ch->tausqinv);
  res->gamma=(ch->gamma);
  res->lambda=(ch->lambda);
  res->pil_pos_mean=(ch->pil_pos_mean);
  res->pil_neg_mean=(ch->pil_neg_mean);
  res->pil_pos_prec=(ch->pil_pos_prec);
  res->pil_neg_prec=(ch->pil_neg_prec);
  res->kap_pos_rate=(ch->kap_pos_rate);
  res->kap_neg_rate=(ch->kap_neg_rate); 
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
  new_vec=(double *) Calloc(len,double);  
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
