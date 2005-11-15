/********************************************/
/*       Authors: Hyungwon Choi             */
/*                Debashis Ghosh            */
/********************************************/
/*           Meta-analysis R pack           */
/********************************************/

#include "R.h"
#include <R_ext/Utils.h>
#include <Rmath.h>
#ifdef WINDOWS
#define fprintf win_print
void win_print(FILE *, char *,...);
#endif
#define MAX_ID 50 /* maximum 50 characters for gene ID */
#define _NA_VAL_ 1e10
#define EPSILON (12*FLT_EPSILON)
#define logit(p) log(p)-log(1-p)

#define _STEPSIZE_ 0.5
#define _SKIP_ 4
#define _BURN_IN_ 500
#define _KAP_MIN_ 3.0
#define _DDFF_ 10.0
#define _PV_GAMMA_ 1.0
#define _POSMIN_ 1e-2


/*********************/
/* gene data struct  */
/*********************/

typedef struct tagARRAY {
  /* char **id;   gene IDs */
  double **d;   /* exprs */
  int nrow;   /* num row */
  int ncol;   /* num col */
  int *label;  /* class label */
  /* char name[MAX_ID]; */ 
} ARRAY;

typedef struct tagARRAY2 {
  /* char **id;   gene IDs */
  double **d;   /* exprs */
  int nrow;   /* num row */
  int ncol;   /* num col */
  int *label;  /* class label */
  double *mean0; /* row sample mean */
  double *mean1;
  double *mean_diff;
  double *var0; /* row sample variance */
  double *var1;
  double *var_sum;
  /* char name[MAX_ID]; */ 
} ARRAY2;

typedef struct tagMARRAY {
  int numdata;
  ARRAY2 *expr;
}  MARRAY;

/***************************************/
/*                  POE                */
/***************************************/
/* structures passed between functions */
/***************************************/

typedef struct tagPR {
  double alpha_mm;
  double alpha_sd;
  double mu_mm;
  double mu_sd;
  double pi_pos_mm;
  double pi_pos_sd;
  double pi_neg_mm;
  double pi_neg_sd;
  double kap_pri_rate;
  double tausqinv_aa;
  double tausqinv_bb;
} PR;

typedef struct tagPP {
  double *alpha_t;
  double *mu_g;
  double *kappa_pos_g;
  double *kappa_neg_g;
  double *sigma_g;
  double *pi_pos_g;
  double *pi_neg_g;
  double mu;
  double tausqinv;
  double gamma;
  double lambda;
  double pil_pos_mean;
  double pil_neg_mean;
  double pil_pos_prec;
  double pil_neg_prec;
  double kap_pos_rate;
  double kap_neg_rate;
  double **poe_mat;
  double **phat_pos;
  double **phat_neg;
  double accept;
} PP;

typedef struct tagCH {
  double **alpha_t;
  double **mu_g;
  double **kappa_pos_g;
  double **kappa_neg_g;
  double **sigma_g;
  double **pi_pos_g;
  double **pi_neg_g;
  double *mu;
  double *tausqinv;
  double *gamma;
  double *lambda;
  double *pil_pos_mean;
  double *pil_neg_mean;
  double *pil_pos_prec;
  double *pil_neg_prec;
  double *kap_pos_rate;
  double *kap_neg_rate;
  double **poe_mat;
  double accept;
} CH;


/*****************************************/
/* This function reads expression matrix */
/* from a double vector and initializes  */
/* all elements to zero                  */
/*****************************************/
void init_ARRAY(double *d, int *nrow, int *ncol, int *label, ARRAY *expr); 
void init_ARRAY2(double *d, int *nrow, int *ncol, int *label, ARRAY2 *expr);
void init_ARRAYS(double *exprs, int *ndata, int *nrow, int *ncol, int *labels, ARRAY2 data[]);
 
/***************************************/
/* This function takes an ARRAY struct */
/* with nrow and ncol initialized at   */
/* init_ARRAY function, dynamically    */
/* allocate the memory space           */
/***************************************/
void malloc_array(ARRAY *expr); 
void malloc_array2(ARRAY2 *expr);

/*****************************************/
/* Following two functions dynamically   */
/* allocates memory for struct PP and CH */
/* only vectors and matrices, no scalars */        
/*****************************************/
void malloc_PP(PP *pp, int *nrow, int *ncol);
void malloc_CH(CH *ch, int *nrow, int *ncol, int *niter);


/*****************************************/
/* Following three functions free memory */
/* dynamically allocated to each struct  */
/*****************************************/
void free_array(ARRAY *expr);
void free_array2(ARRAY2 *expr);
void free_metaarray(MARRAY *data);
void free_PP(PP *pp, int *nrow);
void free_CH(CH *ch, int *nrow, int *ncol, int *niter);


/******************************************/
/* This function converts double array in */
/* matrix form to a vector: matrix is     */
/* considered as an array of rows, e.g.   */
/* so mat[i][j] points to row i and col j */
/******************************************/
void mat2vec(double **mat, double *d, int *nrow, int *ncol);


/*****************************************/
/* Following two functions changes input */
/* from R, which is in vector, into      */
/* proper elements of target structs     */
/*****************************************/
void vec2PP(double *vec, PP *pp, int *nrow, int *ncol);
void vec2PR(double *vec, PR *pr);


/***************************************/
/* Following two functions initializes */
/* all elements of each struct to zero */ 
/***************************************/
void init_PP(PP *pp, int *nrow, int *ncol);
/* void init_CH(CH *ch, int *nrow, int *ncol); */


/*********************************/
/* Binary operations for structs */
/*********************************/
/* void add_PP(PP *pp1, PP *pp2, int *niter); */
/* void divide_PP(PP *pp, int *niter); */


/*****************************************/
/* Converts PP as a numeric vector in    */
/* predetermined fashion, also allocates */
/* memories for PP                       */
/*****************************************/
void PP2vec(double *vec, PP *pp, int *nrow, int *ncol);

/***********************************************/
/* Function need for finalizing report from CH */
/***********************************************/
void update_CH(CH *ch, PP *pp, int iter, int *numiter, int *nrow, int *ncol);
void median_CH(CH *ch, PP *res, int len, int *nrow, int *ncol);
double get_median(double *base, int len);


/***********************************/
/* Functions for Weighted Contrast */
/***********************************/
void weighted_contrast(ARRAY2 data[], int *nd, double *z, int *nrow);
void do_LOWESS(double *x, double *y, int len);
void permute_pval(ARRAY2 data[], int *nd, int *nr, int *nc, int *numperm, double *z, double *p);
void perm(int *in, int *out, int t);
void get_meanvar(ARRAY2 *expr);

/***************************/
/*        Math tools       */
/***************************/
double xlogy(const double x, const double y);
double vec_sum(const double *vec, int len);
double vec_max(const double *vec, int len);
double vec_min(const double *vec, int len);
double vec_mean(const double *vec, int len);
double vec_var(const double *vec, int len);
double log_posterior_gamma(const double a, const double b, const double *sigma_g, const int len);
double log_posterior_kappa(const double kappa, const double *xxx, const int len, const double sigma_g, const double pi_g, const double kr);
void calcor(const double *x, const double *y, const int len, double *num);

/***************************/
/*      Main functions     */
/***************************/
void poe_one_iter(ARRAY *expr, PR *pr, PP *pp);
/* poe_fit(double *expr, int *label, double *prior, double *posterior, int *nrow, int *ncol, int *numiter, double *avgpos); */
/* poe_fit_2(double *expr, int *label, double *prior, double *posterior, int *nrow, int *ncol, int *numiter, double *avgpos); */


