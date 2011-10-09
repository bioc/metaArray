/********************************************/
/*        Authors: Hyungwon Choi            */
/*                 Debashis Ghosh           */
/********************************************/
/*       Probability of Expression          */
/********************************************/
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include "meta.h"

/***********************/
/***********************/
/*   Main functions    */
/***********************/
/***********************/

void poe_fit(double *expr, int *label, double *prior, double *posterior, int *nrow, int *ncol, int *numiter, int *burnin, double *avgpos)
{
  ARRAY data;  
  PR pr;
  PP pp;
  PP res;
  CH ch;
  int k,mm;
  /**********************/
  /*** initialization ***/
  /**********************/
  init_ARRAY(expr,nrow,ncol,label,&data);
  vec2PR(prior, &pr);
  vec2PP(posterior,&pp,nrow,ncol);   
  init_PP(&res,nrow,ncol);
  malloc_CH(&ch,nrow,ncol,numiter);
  Rprintf("%s", "Burn-in\n");
  GetRNGstate();
  /*************************************/
  /*** loop: first burn-in iteration ***/
  /*************************************/
  for(mm=0;mm<*burnin;mm++) {
    poe_one_iter(&data,&pr,&pp);
    if((mm+1)%100==0) Rprintf("%i%s",(mm+1), "\n");
  }
  /**************************************************/
  /*** loop: MC iterate M times while summarizing ***/
  /**************************************************/
  Rprintf("%s", "Main iterations\n");
  for(mm=0;mm<*numiter;mm++) {
    if(_SKIP_>0) {
      for(k=0;k<_SKIP_;k++) poe_one_iter(&data,&pr,&pp);
    }
    poe_one_iter(&data,&pr,&pp);
    if((mm+1)%10 == 0) update_CH(&ch,&pp,mm,numiter,nrow,ncol);    
    if((mm+1)%100 == 0) Rprintf("%i%s",(mm+1), "\n");
  }
  /***********************************************/
  /*** Summarize posterior estimates by median ***/
  /***********************************************/
  Rprintf("%s", "Summary by median\n");
  median_CH(&ch,&res,mm,nrow,ncol);     
  PP2vec(avgpos,&res,nrow,ncol);
  free_array(&data);
  free_PP(&pp,nrow);
  free_PP(&res,nrow);
  free_CH(&ch,nrow,ncol,numiter);
  PutRNGstate();
}


void poe_fit_2(double *expr, int *label, double *prior, double *posterior, int *nrow, int *ncol, int *numiter, double *avgpos)
{
  ARRAY data;  
  PR pr;
  PP pp;
  PP res;
  /* CH ch; */
  int k,mm;
  /**********************/
  /*** initialization ***/
  /**********************/
  init_ARRAY(expr,nrow,ncol,label,&data);
  vec2PR(prior, &pr);
  vec2PP(posterior,&pp,nrow,ncol);   
  init_PP(&res,nrow,ncol);
  /* malloc_CH(&ch,nrow,ncol,numiter); */ 
  GetRNGstate();
  /**************************************************/
  /*** loop: MC iterate M times while summarizing ***/
  /**************************************************/
  /* Rprintf("%s", "Begin Iteration after BurnIn\n"); */
  for(mm=0;mm<*numiter;mm++) {
    if(_SKIP_>0) {
      for(k=0;k<_SKIP_;k++) poe_one_iter(&data,&pr,&pp); 
    }
    poe_one_iter(&data,&pr,&pp);
    /*  update_CH(&ch,&pp,mm,numiter,nrow,ncol);   */ 
    if((mm+1)%100==0) Rprintf("%i%s",(mm+1), " ");
    if((mm+1)%1000==0) Rprintf("%s","\n"); 
  }
  /***********************************************/
  /*** Summarize posterior estimates by median ***/
  /***********************************************/
  PP2vec(avgpos, &pp,nrow,ncol);
  free_array(&data);
  free_PP(&pp,nrow);
  free_PP(&res,nrow);
  PutRNGstate();
}

void poe_one_iter(ARRAY *expr, PR *pr, PP *pp)
{
  /**********************/
  /*** initial values ***/
  /**********************/
  int i, j, ct, tmpint;
  int nr,nc,tt,gg, accept;
  double alpha[expr->ncol];
  double resid_c[expr->ncol];
  double pos_res[expr->ncol];
  double neg_res[expr->ncol];
  double sigma_g[expr->nrow];
  double kappa_pos_g[expr->nrow];
  double kappa_neg_g[expr->nrow];
  double mu_g[expr->nrow];
  double pi_pos_g[expr->nrow];
  double pi_neg_g[expr->nrow];
  double kappa_b, kappa_a, kappa_new, aaa, res_max, shift;         
  double log_prop_old, log_prop_new, log_post_new, log_post_old;
  double d0[expr->nrow][expr->ncol];
  double resid[expr->nrow][expr->ncol];
  double dplus[expr->nrow][expr->ncol];
  double dminus[expr->nrow][expr->ncol];
  double ee[expr->nrow][expr->ncol];
  /* double **phat_pos, **phat_neg; */
  double nnn, post_var, post_mean, post_a, post_b, nn_0;
  double sss[expr->nrow];
  double mmm[expr->nrow];
  double kappa_min, tmp, tmpq;
  double succpos, succneg;
  double tmp1, tmp2, tmp3, tmp4;
  double gamma_new; /*, mu */
  
  /*************************************/
  /*** memory set and initialization ***/
  /*************************************/
  nr=expr->nrow; 
  nc=expr->ncol;
  /* assert(sigma_g=(double *) Calloc(nr,double));
  assert(kappa_pos_g=(double *) Calloc(nr,double));
  assert(kappa_neg_g=(double *) Calloc(nr,double));
  assert(mu_g=(double *) Calloc(nr,double));
  assert(pi_pos_g=(double *) Calloc(nr,double));
  assert(pi_neg_g=(double *) Calloc(nr,double)); */
  for(i=0;i<nc;i++) {
    alpha[i]=pp->alpha_t[i];
  }
  accept=0.0;


  /***************************/
  /***    POE calculation  ***/
  /***     First, memory   ***/
  /***************************/
  /* assert(d0=(double **) Calloc(nr,double *));
  assert(resid=(double **) Calloc(nr,double *));
  assert(dplus=(double **) Calloc(nr,double *));
  assert(dminus=(double **) Calloc(nr,double *)); */
  /*  assert(phat_pos=(double **) Calloc(nr,double *)); */
  /* assert(phat_neg=(double **) Calloc(nr,double *)); */
  /* assert(ee=(double **) Calloc(nr,double *));
  for(i=0;i<nr;i++) {
    assert(d0[i]=(double *) Calloc(nc,double));
    assert(resid[i]=(double *) Calloc(nc,double));
    assert(dplus[i]=(double *) Calloc(nc,double));
    assert(dminus[i]=(double *) Calloc(nc,double)); */
    /*    assert(phat_pos[i]=(double *) Calloc(nc,double)); */
    /*    assert(phat_neg[i]=(double *) Calloc(nc,double)); */
    /* assert(ee[i]=(double *) Calloc(nc,double)); 
  } */

  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      dplus[i][j] = 0.0;
      dminus[i][j] = 0.0;
      d0[i][j] = 0.0;
      ee[i][j] = 0.0;
    }
  }

  for(i=0;i<nr;i++) {
    for(j=0;j<nc;j++) {
      d0[i][j]=(1.0-pp->pi_pos_g[i]-pp->pi_neg_g[i])*(dnorm4(expr->d[i][j],pp->mu_g[i]+pp->alpha_t[j],pp->sigma_g[i],0));      
      resid[i][j]=expr->d[i][j]-pp->mu_g[i]-pp->alpha_t[j];
      if(resid[i][j] >= 0) dplus[i][j]=pp->pi_pos_g[i]*dunif(resid[i][j],0.0,pp->kappa_pos_g[i],0);
      if(resid[i][j] <= 0) dminus[i][j]=pp->pi_neg_g[i]*dunif(resid[i][j],-1.0*pp->kappa_neg_g[i],0.0,0);

      /* Positive */ 
      if(resid[i][j] >= 0) {
        if(dplus[i][j]+d0[i][j] == 0.0) pp->phat_pos[i][j]=0.0;
        else if(dplus[i][j]/(dplus[i][j]+d0[i][j]) > 1.0) pp->phat_pos[i][j]=1.0;
        else if(dplus[i][j]/(dplus[i][j]+d0[i][j]) < 0.0) pp->phat_pos[i][j]=0.0; 
        else pp->phat_pos[i][j]=dplus[i][j]/(dplus[i][j]+d0[i][j]);

        pp->phat_neg[i][j] = 0.0;
      } 
      else {
        if(dminus[i][j]+d0[i][j] == 0.0) pp->phat_neg[i][j]=0.0;
        else if(dminus[i][j]/(dminus[i][j]+d0[i][j]) > 1.0) pp->phat_neg[i][j]=1.0;
        else if(dminus[i][j]/(dminus[i][j]+d0[i][j]) < 0.0) pp->phat_neg[i][j]=0.0; 
        else pp->phat_neg[i][j]=dminus[i][j]/(dminus[i][j]+d0[i][j]); 

        pp->phat_pos[i][j] = 0.0;
      }
      /* ee */
      if(runif(0,1) <= (pp->phat_pos[i][j]+pp->phat_neg[i][j])) ee[i][j] = 1.0;
      else ee[i][j] = 0.0;
      if(resid[i][j] < 0) ee[i][j] = (-1.0) * ee[i][j] ;
    }
  }
  for(j=0;j<nc;j++) {
    if(expr->label[j]==1) {
      for(i=0;i<nr;i++) {
        ee[i][j]=0.0;
      }  	  
    }
  }

  /************************/
  /***  alpha and mu    ***/
  /************************/
  /***      alpha       ***/
  /************************/

  /****************/
  /***     mu   ***/
  /****************/
  for(gg=0;gg<nr;gg++) {
    nnn=0.0;
    for(j=0;j<nc;j++) nnn +=(1.0 - fabs(ee[gg][j]));

    if(!ISNA(nnn) && !ISNAN(nnn) && R_FINITE(nnn) && (nnn>0.0)) {
      ct=0;
      for(j=0;j<nc;j++) ct+=((int) (ee[gg][j]==0.0));
      /* mmm=(double *) Calloc(ct,double); */
      ct=0;
      for(j=0;j<nc;j++) {
        if(ee[gg][j]==0.0) {
          mmm[ct]=(expr->d[gg][j]-alpha[j])/pow(pp->sigma_g[gg],2);
	  ct++;
	}
      }
      tmp1=0.0;
      for(i=0;i<ct;i++) {
        tmp1+=mmm[i];
      }
      tmp1=tmp1/((double) ct);
      post_var=1.0/(nnn/pow(pp->sigma_g[gg],2.0)+1.0/pp->tausqinv);
      post_mean=(nnn*tmp1+pp->mu/pow(pp->tausqinv,2.0))*post_var;
      mu_g[gg] = rnorm(post_mean,sqrt(post_var));
      /* mu_g[gg] = ((!R_FINITE(tmp1)) || (ISNA(tmp1)) || (ISNAN(tmp1)) ? pp->mu_g[gg] : tmp1); */
      /* Free(mmm); */ 
    }
    else {
      mu_g[gg]=pp->mu_g[gg]; 
    }
    nn_0=0.0;
    for(j=0;j<nc;j++) nn_0+=(1.0-fabs(ee[gg][j]));
    post_a=pp->gamma+nn_0/2.0;
    tmp3=0.0;
    tmp1=0.0;
    for(j=0;j<nc;j++) {
      tmp1=1.0-fabs(ee[gg][j]);
      tmp2=expr->d[gg][j]-mu_g[gg]-alpha[j];
      tmp3+=tmp1*tmp2*tmp2;
    }
    post_b=pp->lambda+0.5*tmp3;
    kappa_min=fmin2(kappa_pos_g[gg],kappa_neg_g[gg]);
    tmp1=pgamma(pow(_KAP_MIN_/kappa_min,2.0),post_a,1.0/post_b,1.0,0);
    tmp=runif(tmp1,1.0);
    tmpq=qgamma(tmp,post_a,1.0/post_b,1.0,0);
    sigma_g[gg]= ((R_FINITE(tmpq) && (!ISNA(tmpq)) && (!ISNAN(tmpq))) ? sqrt(1.0/tmpq) : pp->sigma_g[gg]); 
  }    
    /* if(ISNA(sigma_g[gg]) || ISNAN(sigma_g[gg])) 
       print error to stderr */

  for(tt=0;tt<nc;tt++) {
    nnn=0.0;

    for(i=0;i<nr;i++) nnn += (1.0 - fabs(ee[i][tt]));
    if((!ISNA(nnn)) && (!ISNAN(nnn)) && (R_FINITE(nnn)) && (nnn > 0.0)) {
      ct=0;
      for(i=0;i<nr;i++) ct+=((int) (ee[i][tt]==0.0));
      /* sss=(double *) Calloc(ct,double);
	 mmm=(double *) Calloc(ct,double); */
      ct=0;
      for(i=0;i<nr;i++) {
        if(ee[i][tt]==0.0) {
          tmp2=pow(pp->sigma_g[i],2.0);
          sss[ct]=1.0/tmp2;
          mmm[ct]=(expr->d[i][tt]-pp->mu_g[i])/tmp2;
          ct++;
        }
      }
      tmp1=vec_sum(sss,ct);
      tmp2=vec_mean(mmm,ct);
      post_var = 1.0/(tmp1+(1.0/pow(pr->alpha_sd,2.0)));
      post_mean = (nnn*tmp2 + pr->alpha_mm/pow(pr->alpha_sd,2.0)) * post_var;
      tmp1 = rnorm(post_mean,sqrt(post_var));
      alpha[tt] = ((!R_FINITE(tmp1)) || (ISNA(tmp1)) || (ISNAN(tmp1)) ? pp->alpha_t[tt] : tmp1);
      /* Free(sss);
      Free(mmm); */
    }
    else {
      alpha[tt]=pp->alpha_t[tt];
    }
  }
  tmp2 = 0.0;
  for(i=0;i<nc;i++) tmp2 += alpha[i];
  tmp2 = tmp2 / ((double) nc);
  for(i=0;i<nc;i++) alpha[i] -= tmp2;

  /**********************/  
  /*** Sample Kappa's ***/
  /**********************/
  for(gg=0;gg<nr;gg++) 
  {
    for(j=0;j<nc;j++) resid_c[j] = (expr->d[gg][j])-(pp->alpha_t[j])-(pp->mu_g[gg]);
    /*********************/
    /*** Positive Side ***/
    /*********************/
    kappa_b=pp->kap_pos_rate; 
    kappa_a=1.0;
    ct=0;
    for(j=0;j<nc;j++) {
      if(resid_c[j]>0) ct++;
    }
    /* assert(pos_res=(double *) Calloc(ct,double)); */
    ct=0;
    res_max = 0.0;
    for(j=0;j<nc;j++) {
      if(resid_c[j]>=0) {
        pos_res[ct]=resid_c[j];
        if(resid_c[j] >= res_max) res_max = resid_c[j]*1.01;
        ct++;
      }
    }    
    
    if(res_max > _KAP_MIN_ * pp->sigma_g[gg]) { 
      /* kappa_new= res_max + rgamma(kappa_a,1.0/kappa_b); */
      kappa_a = 1.0 + pp->kap_pos_rate * (res_max - _KAP_MIN_ * pp->sigma_g[gg]);
    }
    shift = fmax2(res_max, _KAP_MIN_ * pp->sigma_g[gg]);
    kappa_new = res_max + rgamma(kappa_a, 1.0/kappa_b);
    log_prop_old=dgamma(pp->kappa_pos_g[gg],kappa_a,1.0/kappa_b,0);
    log_prop_new=dgamma(kappa_new,kappa_a,1.0/kappa_b,0);
    log_post_new=log_posterior_kappa(kappa_new,pos_res,ct,pp->sigma_g[gg],pp->pi_pos_g[gg],pp->kap_pos_rate);
    log_post_old=log_posterior_kappa(pp->kappa_pos_g[gg],pos_res,ct,pp->sigma_g[gg],pp->pi_pos_g[gg],pp->kap_pos_rate);
    if(ISNAN(log_post_old) || ISNAN(log_post_old)) aaa=0.0;
    if(ISNAN(log_post_new) || ISNAN(log_post_new)) aaa=0.0;
    aaa=exp(log_post_new-log_post_old-log_prop_new+log_prop_old);
    if(ISNA(aaa)||ISNAN(aaa)||!R_FINITE(aaa)) aaa=0.0;
    if(aaa > runif(0,1) && !ISNAN(kappa_new) && !ISNA(kappa_new) && R_FINITE(kappa_new)) { 
      kappa_pos_g[gg]=kappa_new;
    }
    else {
      kappa_pos_g[gg]=pp->kappa_pos_g[gg];
    }
    /* Free(pos_res); */ 

    /*********************/
    /*** Negative side ***/
    /*********************/
    kappa_b=pp->kap_neg_rate;
    kappa_a=1.0;
    ct=0;
    for(j=0;j<nc;j++) if(resid_c[j]<0) ct++;
    /* assert(neg_res=(double *) Calloc(ct,double)); */
    ct=0;
    res_max = 0.0;
    for(j=0;j<nc;j++) {
      if(resid_c[j]<=0) {
        neg_res[ct]=-1.0*resid_c[j];
        if(resid_c[j] >= res_max) res_max = -1.0*resid_c[j]*1.01;
        ct++;
      }
    }    
    if(res_max > _KAP_MIN_*pp->sigma_g[gg]) { 
      kappa_a = 1.0 + pp->kap_neg_rate * (res_max - _KAP_MIN_ * pp->sigma_g[gg]);
    }
    shift = fmax2(res_max, _KAP_MIN_ * pp->sigma_g[gg]);
    kappa_new = res_max + rgamma(kappa_a, 1.0/kappa_b);
    log_prop_old=dgamma(pp->kappa_neg_g[gg],kappa_a,1.0/kappa_b,0);
    log_prop_new=dgamma(kappa_new,kappa_a,1.0/kappa_b,0);
    log_post_new=log_posterior_kappa(kappa_new,neg_res,ct,pp->sigma_g[gg],pp->pi_neg_g[gg], pp->kap_neg_rate);
    log_post_old=log_posterior_kappa(pp->kappa_neg_g[gg],neg_res,ct,pp->sigma_g[gg],pp->pi_neg_g[gg],pp->kap_neg_rate);
    if(ISNA(log_post_old) || ISNAN(log_post_old)) aaa=0.0;
    if(ISNAN(log_post_new) || ISNAN(log_post_new)) aaa = 0.0;
    aaa=exp(log_post_new-log_post_old-log_prop_new+log_prop_old);
    if(ISNA(aaa)||ISNAN(aaa)||!R_FINITE(aaa)) aaa=0.0;
    if(aaa>runif(0,1) && !ISNA(kappa_new) && !ISNAN(kappa_new) && R_FINITE(kappa_new)) {
       kappa_neg_g[gg]=kappa_new;
    }
    else {
       kappa_neg_g[gg]=pp->kappa_neg_g[gg];
    }
    /* Free(neg_res); */
  }  


  for(gg=0;gg<nr;gg++) {
    succpos=0.0;
    succneg=0.0;
    for(j=0;j<nc;j++) {
      succpos +=  (ee[gg][j] > 0.0);
      succneg +=  (ee[gg][j] < 0.0);
    }

    pi_pos_g[gg]=rbeta(succpos+1.0, ((double) nc)-succpos+1.0);
    pi_neg_g[gg]=rbeta(succneg+1.0, ((double) nc)-succneg+1.0);

    tmp1=dnorm4(logit(pi_pos_g[gg]),pp->pil_pos_mean,sqrt(pp->pil_pos_prec),0);
    tmp2=dnorm4(logit(pp->pi_pos_g[gg]),pp->pil_pos_mean,sqrt(pp->pil_pos_prec),0);
    aaa=tmp1/tmp2;
    if(ISNAN(aaa) || ISNA(aaa) || !R_FINITE(aaa)) aaa=0.0;
    tmpint=0;
    for(j=0;j<nc;j++) tmpint += expr->label[j];
    if(tmpint==0 && pi_pos_g[gg]>0.5) pi_pos_g[gg]=pp->pi_pos_g[gg];
    if(aaa<runif(0,1)) pi_pos_g[gg]=pp->pi_pos_g[gg];    
    tmp1=dnorm4(logit(pi_neg_g[gg]),pp->pil_neg_mean,sqrt(pp->pil_neg_prec),0);
    tmp2=dnorm4(logit(pp->pi_neg_g[gg]),pp->pil_neg_mean,sqrt(pp->pil_neg_prec),0);
    aaa=tmp1/tmp2;
    if(ISNAN(aaa) || ISNA(aaa) || !R_FINITE(aaa)) aaa=0.0;
    if(tmpint==0 && pi_neg_g[gg] > 0.5) pi_neg_g[gg]=pp->pi_neg_g[gg];
    if(aaa<runif(0,1) || pi_pos_g[gg]+pi_neg_g[gg] > 1.0) pi_neg_g[gg]=pp->pi_neg_g[gg];    
  } 


  /******************************/
  /*** Hyperparameters sample ***/
  /******************************/
  post_var=1.0/(((double) nr)*pp->tausqinv+1.0/pow(pr->mu_sd,2));
  post_mean=(((double) nr)*pp->tausqinv*(vec_mean(mu_g,nr))+(pr->mu_mm/pow(pr->mu_sd,2)))*post_var;
  pp->mu=rnorm(post_mean,sqrt(post_var));
  post_a=pr->tausqinv_aa+(((double) nr)/2.0);
  tmp1=0.0;
  for(i=0;i<nr;i++) tmp1+=pow(mu_g[i]-pp->mu,2.0);
  post_b=pr->tausqinv_bb+0.5*tmp1;
  pp->tausqinv=rgamma(post_a,1.0)/post_b;
  /*pp->tausqinv = (!R_FINITE(tmp) || ISNA(tmp) || ISNAN(tmp) ? pp->tausqinv : tmp);*/
  pp->kap_pos_rate=rgamma(((double)(nr+1)),1.0/(pr->kap_pri_rate+vec_sum(kappa_pos_g,nr)));
  /* pp->kap_pos_rate=(!R_FINITE(tmp) || ISNA(tmp) || ISNAN(tmp) ? pp->kap_pos_rate : tmp); */
  pp->kap_neg_rate=rgamma(((double)(nr+1)),1.0/(pr->kap_pri_rate+vec_sum(kappa_neg_g,nr)));
  /* pp->kap_neg_rate=(!R_FINITE(tmp) || ISNA(tmp) || ISNAN(tmp) ? pp->kap_neg_rate : tmp); */

  tmp1=0.0;
  tmp2=0.0;
  tmp3=0.0;
  tmp4=0.0;
  for(i=0;i<nr;i++) {
    tmp1+=logit(pi_pos_g[i]);
    tmp2+=logit(pi_neg_g[i]);
    tmp3+=pow(logit(pi_pos_g[i])-pp->pil_pos_mean,2.0);
    tmp4+=pow(logit(pi_neg_g[i])-pp->pil_neg_mean,2.0);
  } 
  tmp1=tmp1/((double) nr);
  tmp2=tmp2/((double) nr);

  post_var=1.0/(((double) nr)*pp->pil_pos_prec+1.0/pow(pr->mu_sd,2.0));
  post_mean=(((double) nr)*pp->pil_pos_prec*tmp1)+(pr->pi_pos_mm/pow(pr->pi_pos_sd,2.0))*post_var;
  pp->pil_pos_mean = rnorm(post_mean,sqrt(post_var));
  /* pp->pil_pos_mean = (!R_FINITE(tmp) || ISNA(tmp) || ISNAN(tmp) ? pp->pil_pos_mean : tmp) ; */

  post_a=pr->tausqinv_aa+(((double) nr) / 2.0);
  post_b=pr->tausqinv_bb+0.5*tmp3;
  pp->pil_pos_prec = rgamma(post_a,1.0)/post_b;
  /* pp->pil_pos_prec = (!R_FINITE(tmp) || ISNA(tmp) || ISNAN(tmp) ? pp->pil_pos_prec : tmp); */

  post_var=1.0/(((double) nr)*pp->pil_neg_prec+1.0/pow(pr->mu_sd,2.0));
  post_mean=(((double) nr)*pp->pil_neg_prec*tmp2)+(pr->pi_neg_mm/pow(pr->pi_neg_sd,2.0))*post_var;
  pp->pil_neg_mean = rnorm(post_mean,sqrt(post_var));
  /* pp->pil_neg_mean = (!R_FINITE(tmp) || ISNA(tmp) || ISNAN(tmp) ? pp->pil_neg_mean : tmp); */

  post_a=pr->tausqinv_aa+(((double) nr) / 2.0);
  post_b=pr->tausqinv_bb+0.5*tmp4;
  pp->pil_neg_prec = rgamma(post_a,1.0)/post_b;
  /* pp->pil_neg_prec = (!R_FINITE(tmp) || ISNA(tmp) | ISNAN(tmp) ? pp->pil_neg_prec : tmp); */

  tmp1=0.0;
  for(i=0;i<nr;i++) tmp1+=1.0/pow(sigma_g[i],2.0);
  tmp=rgamma(((double) nr)*pp->gamma+1.0,1.0/tmp1);
  pp->lambda = (!R_FINITE(tmp) || ISNA(tmp) || ISNAN(tmp) ? pp->lambda : tmp);
  gamma_new=pp->gamma+_STEPSIZE_*sqrt(_PV_GAMMA_)*rt(_DDFF_);     

  log_post_new=log_posterior_gamma(gamma_new,pp->lambda,sigma_g,nr);
  log_post_old=log_posterior_gamma(pp->gamma,pp->lambda,sigma_g,nr);

  aaa=exp(log_post_new-log_post_old);

  if(ISNAN(log_post_old) || ISNAN(log_post_old) || !R_FINITE(log_post_old)) aaa=0.0;
  if(ISNAN(log_post_new) || ISNAN(log_post_new) || !R_FINITE(log_post_new)) aaa=0.0;
  if(ISNAN(aaa) || ISNA(aaa) || !R_FINITE(aaa)) aaa=0.0;
  if(aaa>runif(0,1)) {
    pp->gamma=gamma_new;
    pp->accept=1.0;
  }
  else {
    pp->gamma=pp->gamma;
    pp->accept=0.0;
  }

  /***************************/
  /***  Update POE matrix  ***/
  /***************************/
  for(j=0;j<nc;j++) pp->alpha_t[j]=alpha[j];
  for(i=0;i<nr;i++) {
    pp->mu_g[i]=mu_g[i];
    pp->sigma_g[i]=sigma_g[i];
    pp->pi_pos_g[i]=pi_pos_g[i];
    pp->pi_neg_g[i]=pi_neg_g[i];
    pp->kappa_pos_g[i]=kappa_pos_g[i];
    pp->kappa_neg_g[i]=kappa_neg_g[i];
    for(j=0;j<nc;j++) {
      pp->poe_mat[i][j]=pp->phat_pos[i][j]-pp->phat_neg[i][j];
      /* if(ISNAN(pp->poe_mat[i][j]) || ISNA(pp->poe_mat[i][j]) || !R_FINITE(pp->poe_mat[i][j])) Rprintf("NA in POE\n"); */
      /*  Rprintf("%f\t", pp->poe_mat[i][j]); */
      /* pp->phat_pos[i][j]=phat_pos[i][j];  
	 pp->phat_neg[i][j]=phat_neg[i][j];  */
    }
  } 

  /**************************/
  /*** Garbage collection ***/
  /**************************/
  /* Free(resid_c);
  Free(pi_neg_g);
  Free(pi_pos_g);
  Free(mu_g);
  Free(kappa_neg_g);
  Free(kappa_pos_g);
  Free(sigma_g);
  Free(alpha);
  for(i=0;i<nr;i++) {
    Free(d0[i]);
    Free(resid[i]);
    Free(dplus[i]);
    Free(dminus[i]);
    Free(ee[i]);
  }
  Free(d0); 
  Free(resid);
  Free(dplus);
  Free(dminus);
  Free(ee); */
}




