// SPDX-License-Identifier: BSD-3-Clause
/* optim.h — 数值优化模块头文件 */
#ifndef OPTIM_H_
#define OPTIM_H_

#include "types.h"
#include "linalg.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --- lnsrchmp.h --- */
#define EPSILON 2.7182818284590452353602874713526624977572

int stopcheck_mt(double fx, int N, double *xc, double *xf, double *jac, double *dx, double fsval, double gtol, double stol, int retval);

int stopcheck2_mt(double fx, int N, double fo, double *jac, double *dx, double eps,double stoptol, double functol, int retval);

int stopcheck3_mt(double *xi,double *xf,double fx, int N, double fo, double *jac, double *dx, double eps,
		double stoptol, double functol, int retval);

int grad_fd(custom_function *funcpt,custom_gradient *funcgrad, double *x, int N, double *dx, double eps2, double *f);

int grad_cd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
		double eps3, double *f);

int grad_calc2(custom_function *funcpt, double *x, int N, double *dx, double eps3, double *f);

int grad_calc(custom_function *funcpt, double *x, int N, double *dx, double eps2, double *f);

int cstep(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy, double *stp, double *fp, double *dp, int *brackt,
	double  stpmin, double stpmax);

int cvsrch(custom_function *funcpt, custom_gradient *funcgrad, double *x, double *f, double *g, double *stp, double *s, int N, double *dx, double maxstep,
	int MAXITER,double eps2,double ftol, double gtol, double xtol);

int lnsrchmt(custom_function *funcpt, custom_gradient *funcgrad, double *xi, double *f, double *jac, double *alpha, double *p, int N, double *dx, double maxstep, int MAXITER,
		double eps2,double ftol, double gtol, double xtol, double *x);

/* --- newtonmin.h --- */
double signx(double x);

double cholmod(double *A, int N, double *L, double eps,double maxinp);

double modelhess(double *A,int N,double *dx,double eps,double *L);

void linsolve_lower(double *L,int N,double *b,double *x);

int hessian_fd(custom_function *funcpt, double *x, int N, double *dx, double eps, double *f);

int hessian_fd2(custom_function *funcpt,double *x,int N,double *dx,double eps,double *f);

void fdjac(custom_gradient *funcgrad, double *x, int N, double *jac, double *dx, double eps2, double *J);

void hessian_fdg(custom_gradient *funcgrad, double *x, int N, double *jac, double *dx, double eps2, double *H);

int hessian_opt(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *jac,
		double *dx,double eps,double eps2,double *H);

int lnsrch(custom_function *funcpt, double *xi, double *jac, double *p, int N, double * dx, double maxstep, double stol, double *x);

int lnsrchmod(custom_function *funcpt, custom_gradient *funcgrad, double *xi, double *jac, double *p, int N, double * dx, double maxstep,
	double eps2,double stol,double *x,double *jacf);
	
int lnsrchcg(custom_function *funcpt, custom_gradient *funcgrad, double *xi, double *jac, double *p, int N, double * dx, double maxstep,
	double eps2,double stol,double *x,double *jacf);
	
int stopcheck(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol,int retval);

int newton_min_func(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval, double maxstep, int MAXITER,
	int *niter, double eps, double gtol, double stol, double *xf);

int trsrch(custom_function *funcpt, double *xi, double *jac, double *sN, int N, double * dx, double maxstep,
		int iter,double *L,double *hess,double stol,double *ioval,double eps,double *x);
		
void trstep(double *jac,double *sN,int N,double * dx,double *L,double *hess,double nlen,double *ioval,double eps,
		double *step);

int trupdate(custom_function *funcpt, double *xi, double *jac, double *step, int N, double * dx, double maxstep,
		int retcode,double *L,double *hess,double stol,int method,double *ioval,double *xprev,double *funcprev,double *x);		

void trstep_ddl(double *jac,double *sN,int N,double * dx,double maxstep,double *L,double *hess,double nlen,double *ioval,
		double *ssd,double *v,double *step);		
		
int trsrch_ddl(custom_function *funcpt, double *xi, double *jac, double *sN, int N, double * dx, double maxstep,
		int iter,double *L,double *hess,double stol,double *ioval,double *x);		

int newton_min_trust(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval, double delta,
		int method,int MAXITER,int *niter,double eps,double gtol,double stol,double *xf);

/* --- conjgrad.h --- */
int ichol(double *A, int N);

int stopcheck2(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol) ;

int cgpr_mt(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx,double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf); //Polak Ribiere + (More Thuentes Line Search)

int conjgrad_min_lin(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf);

/* --- secant.h --- */
void bfgs_naive(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf);

int bfgs_min_naive(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval,double maxstep, int MAXITER,
		double eps, double *xf);

void bfgs_factored(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf);

int bfgs_min(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval, double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double stol,double *xf);

int bfgs_min2(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, int m, double *dx, double fsval, double maxstep, int MAXITER, int *niter,
	double eps, double gtol, double ftol, double xtol, double *xf);

void inithess_l(double *H, int N, int k, double *tsk, double *tyk, double *dx);

void bfgs_rec(double *H, int N, int iter, int m, double *jac, double *sk, double *yk, double *r);

int bfgs_l_min(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, int m, double *dx, double fsval,double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf);

/* --- brent.h --- */
// see https://people.math.sc.edu/Burkardt/f_src/brent/brent.html
double brent_zero(custom_funcuni *funcuni,double a, double b, double tol, double eps);

double brent_local_min(custom_funcuni *funcuni, double a, double b, double t, double eps, double *x);

/* --- neldermead.h --- */
/*
 * neldormead.h
 *
 *  Created on: Jan 5, 2014
 *      Author: HOME
 */
int nel_min(custom_function *funcpt,double *xc,int N,double *dx,double fsval,int MAXITER,int *niter,
		double eps,double *xf);

/* --- optimc.h --- */
/*
 * optimc.h
 *
 *  Created on: Mar 16, 2014
 *      Author: HOME
 */
typedef struct opt_set* opt_object;

opt_object opt_init(int N);

struct opt_set{
	int N;
	double objfunc;
	double eps;
	double gtol;
	double stol;
	double ftol;
	double xtol;
	double maxstep;
	int MaxIter;
	int Iter;
	int Method;
	int retval;
	char MethodName[50];
	double xopt[1];
};

typedef struct nls_set* nls_object;

nls_object nls_init(int M,int N);

struct nls_set{
	int M;
	int N;
	double eps;
	double epsfcn;
	double factor;
	double gtol;
	double ftol;
	double xtol;
	int MaxIter;
	int Maxfev;
	int Iter;
	int nfev;
	int njev;
	int ldfjac;
	int mode;
	int retval;
	double xopt[1];
};

void setnlsTOL(nls_object obj,double gtol,double ftol,double xtol);

void optsummary(opt_object obj);

void setMaxIter(opt_object obj,int MaxIter);

void setMaxStep(opt_object obj, double maxstep);

void setTOL(opt_object obj,double gtol,double stol,double ftol,double xtol);

int fminsearch(custom_function *funcpt,int N,double *xi,double *xf);

double fminbnd(custom_funcuni *funcuni,double a, double b);

int fminunc(custom_function *funcpt,custom_gradient *funcgrad,int N,double *xi,double maxstep, int method,double *xf);

int fminnewt(custom_function *funcpt, custom_gradient *funcgrad, int N, double *xi,
	     double delta,double *dx,double fsval,double maxstep,int method,double *xf);

double brent_local_min(custom_funcuni *funcuni,double a, double b, double t, double eps, double *x);

void optimize(opt_object obj, custom_function *funcpt, custom_gradient *funcgrad, int N, double *xi,
		int method);

void free_opt(opt_object object);

int levmar(custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *xi,int M, int N,double *xf);

void nls(nls_object obj, custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *xi);

void nls_scale(nls_object obj, custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *diag,double *xi);

void nlssummary(nls_object obj);

void free_nls(nls_object object);

#ifdef __cplusplus
}
#endif

#endif /* OPTIM_H_ */

