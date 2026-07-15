// SPDX-License-Identifier: BSD-3-Clause
/* linalg.h — 矩阵与最小二乘模块头文件 */
#ifndef LINALG_H_
#define LINALG_H_

#include "cstats.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --- matrix.h --- */
/*
 * matrix.h
 *
 *  Created on: Jul 1, 2013
 *      Author: USER
 */

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define CUTOFF 192
#define TOL 1e-12
#define BLOCKSIZE 64
#define TBLOCK 64
#define SVDMAXITER 50
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define RSVD_POWER_ITERATIONS 5


double macheps();

double pmax(double a, double b);

double pmin(double a, double b);

int imax(int a, int b);

int imin(int a, int b);

double signx(double x);

int iabs(int a);

double l2norm(double *vec, int N);

int compare (const void* ind1, const void* ind2);

int compare_ascending (const void* ind1, const void* ind2);

void sort1d(double* v,int N, int* pos);

void sort1d_ascending(double* v,int N, int* pos);

//Array Parallel Implementation may have a lot of overhead

double array_max_abs(double *array,int N);

double array_max(double *array,int N);

double array_min(double *array,int N);

//void mmult(double* A, double *B, double *C,int ra,int ca, int rb, int cb);

void dtranspose(double *sig, int rows, int cols,double *col);

void stranspose(double *sig, int rows, int cols,double *col);

void rtranspose(double *m, int rows, int cols,double *n, int r, int c);

void ctranspose(double *sig, int rows, int cols,double *col);

void mtranspose(double *sig, int rows, int cols,double *col);

void itranspose(double *A, int M, int N);

//int minverse(double *xxt, int p);

void mdisplay(double *A, int row, int col);

void madd(double* A, double* B, double* C,int rows,int cols);

void msub(double* A, double* B, double* C,int rows,int cols);

void scale(double *A, int rows, int cols, double alpha);

void nmult(double* A, double* B, double* C,int m,int n, int p);

void tmult(double* A, double* B, double* C,int m,int n, int p);

void recmult(double* A, double* B, double* C,int m,int n, int p,int sA,int sB, int sC);

void rmult(double* A, double* B, double* C,int m,int n, int p);

int findrec(int *a, int *b, int *c);

double house_2(double*x,int N,double *v);

void add_zero_pad(double *X, int rows, int cols, int zrow, int zcol,double *Y);

void remove_zero_pad(double *X, int rows, int cols, int zrow, int zcol,double *Y);

void madd_stride(double* A, double* B, double* C,int rows,int cols,int sA,int sB,int sC);

void msub_stride(double* A, double* B, double* C,int rows,int cols,int sA,int sB,int sC);

void rmadd_stride(double* A, double* B, double* C,int rows,int cols,int p,int sA,int sB,int sC);

void rmsub_stride(double* A, double* B, double* C,int rows,int cols,int p,int sA,int sB,int sC);

void srecmult(double* A, double* B, double* C,int m,int n, int p,int sA,int sB,int sC);

void smult(double* A, double* B, double* C,int m,int n, int p);

void mmult(double* A, double* B, double* C,int m,int n, int p);

void ludecomp(double *A,int N,int *ipiv);

int rludecomp(double *A, int M, int N, int *ipiv);

void getPLU(double *A, int M , int N, int *ipiv,double *P, double *L, double *U);

void getPU(double *A, int M, int N, int *ipiv, double *P,double *U);

double* marsaglia_generate(double *values, int N, double average, double  deviation);

void random_matrix(double *A, int M, int N);

void linsolve(double *A,int N,double *b,int *ipiv,double *x);

void minverse(double *A,int M,int *ipiv,double *inv);

void eye(double *mat,int N);

void eye_scale(double *mat, int N, double lambda);

double house(double*x,int N,double *v);

void housemat(double *v, int N,double beta,double *mat);

void qrdecomp(double *A, int M, int N,double *bvec);

void getQR(double *A,int M,int N,double *bvec,double *Q, double *R);

void hessenberg(double *A,int N);

void francisQR(double *A,int N);

void eig22(double *A, int stride,double *eigre,double *eigim);

int francis_iter(double *A, int N, double *H);

void eig(double *A,int N,double *eigre,double *eigim);

void eigensystem(double *mat, int N, double *eval, double *evec);;

int cholu(double *A, int N);

int bcholu(double *A, int N);

int chol(double *A, int N);

void chold(double *A, int N);

void svd_sort(double *U,int M,int N,double *V,double *q);

int svd(double *A,int M,int N,double *U,double *V,double *q);

int svd_transpose(double *A, int M, int N, double *U, double *V, double *q);

int rank(double *A, int M,int N);

void rsvd(double *A, int M, int N,int K, int oversample, int n_iter,double *U, double *V, double *S);

/* --- lls.h --- */
/*
 * lls.h
 *
 *  Created on: Apr 14, 2014
 *      Author: HOME
 */
int lls_normal(double *A,double *b,int M,int N,double *x);

int lls_qr(double *A,double *b,int M,int N,double *xo);

void bidiag(double *A, int M, int N);

void bidiag_orth(double *A, int M, int N,double *U,double *V);

int svd_gr(double *A,int M,int N,double *U,double *V,double *q);

int svd_gr2(double *A,int M,int N,double *U,double *V,double *q);

int minfit(double *AB,int M,int N,int P,double *q);

int lls_svd(double *Ai,double *bi,int M,int N,double *xo);

int lls_svd2(double *Ai,double *bi,int M,int N,double *xo);

/* --- regression.h --- */
/*
 * regression.h
 *
 *  Created on: Jun 5, 2013
 *      Author: USER
 */
typedef struct bparam_t {
  double value;
  double lower;
  double upper;
  double stdErr;
} bparam;

typedef struct reg_set* reg_object;

reg_object reg_init(int N, int p);

struct reg_set{
	int N;
	int p;
	double alpha;
	double sigma;
	double sigma_lower;
	double sigma_upper;
	double r2;
	double r2adj;
	double R2[2];
	char lls[10];
	int df;
	int intercept;
	int rank; // Use method "qr" to calculate the rank
	double TSS;
	double ESS;
	double RSS;
	int df_ESS;
	int df_RSS;
	double FStat;
	double PVal;
	double loglik;
	double aic;
	double bic;
	double aicc;
	bparam beta[1];
};

void linreg_clrm(double *x,double *y, int N, double* b,
		double *var,double *res,double alpha,double *anv,
		double* ci_lower, double* ci_upper);
		
void zerohyp_clrm(int N,double *b, double *val, double *tval, double *pval);

//void linreg_multi2(int p, double *x,double *y, int N, double* b); // p number of variables.
// p = 2 for one dependent variable	and one independent variable
// p = 3 for one dependent variable	and two independent variables etc.

void linreg_multi(int p, double *xi,double *y, int N, double* b,double *sigma2,
			double *xxti,double *R2,double *res,double alpha,double *anv,
			double* ci_lower, double* ci_upper,int *rank,char *llsmethod, int intercept);
		
void zerohyp_multi(int N,double *b,int p, double *varcovar, double *tval, double *pval);

void regress(reg_object obj,double *x,double *y,double *res,double *varcovar,double alpha);

void regress_poly(reg_object obj,double *x,double *y,double *res,double *varcovar,double alpha);

void setIntercept(reg_object obj,int intercept);

void setLLSMethod(reg_object obj,char *llsmethod);

void summary(reg_object obj);

void anova(reg_object obj);

void anova_list(reg_object *list, int N);

void confint(reg_object obj);

void zerohyp_val(reg_object obj, double *tval, double *pval);

double fitted(reg_object obj,double *inp,double *varcovar,double *var);

void free_reg(reg_object object);

/* --- nls.h --- */
/*
 * nls.h
 *
 *  Created on: May 21, 2014
 *      Author: HOME
 */
double enorm(double *x, int N);

void qrfac(double *A, int M, int N, int lda, int pivot, int *ipvt, int lipvt,double *rdiag, double *acnorm,double eps);

void qrsolv(double *r,int ldr,int N,int *ipvt,double *diag,double *qtb,double *x,double *sdiag);

void fdjac2(custom_funcmult *funcmult, double *x, int M, int N, double *fvec, double *fjac, int ldfjac,
		double epsfcn,double eps);

void lmpar(double *r,int ldr,int N,int *ipvt,double *diag,double *qtb,double delta,double *par,double *x,double *sdiag);

int lmder(custom_funcmult *funcmult,custom_jacobian *jacobian,double *xi,int M, int N,
		double *fvec,double *fjac,int ldfjac,int maxfev,double *diag,int mode,double factor,int nprint,
		double eps,double ftol,double gtol,double xtol,int *nfev,int *njev,int *ipvt, double *qtf);

int lmdif(custom_funcmult *funcmult, double *x, int M, int N, double *fvec, double *fjac, int ldfjac,
		int maxfev,double *diag,int mode,double factor,int nprint,double eps,double epsfcn,double ftol,double gtol,
		double xtol,int *nfev,int *njev,int *ipvt, double *qtf);

#ifdef __cplusplus
}
#endif

#endif /* LINALG_H_ */

