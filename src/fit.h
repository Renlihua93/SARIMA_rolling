// SPDX-License-Identifier: BSD-3-Clause
/*
 * fit.h — 季节 ARIMA 拟合与预测（精简）
 */
#ifndef FIT_H_
#define FIT_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct alik_css_seas_set* alik_css_seas_object;

alik_css_seas_object alik_css_seas_init(int p, int d, int q, int s, int P, int D, int Q, int N);

struct alik_css_seas_set{
	int p;
	int d;
	int q;
	int s;
	int P;
	int D;
	int Q;
	int r;
	int pq;
	int length;
	int N;
	int M;
	double eps;
	double mean;
	double ssq;
	double loglik;
	int offset;
	double x[1];
};

typedef struct alik_seas_set* alik_seas_object;

alik_seas_object alik_seas_init(int p, int d, int q, int s, int P, int D, int Q, int N);

struct alik_seas_set{
	int p;
	int d;
	int q;
	int s;
	int P;
	int D;
	int Q;
	int r;
	int pq;
	int length;
	int N;
	int M;
	double eps;
	double mean;
	double ssq;
	double loglik;
	int offset;
	double x[1];
};

int starma(int ip, int iq, double *phi, double *theta, double *A, double *P, double *V);

void karma(int ip, int iq, double *phi, double *theta, double *A, double *P, double*V, int N,
	double *W, double *resid, double *sumlog, double *ssq, int iupd, double delta, int *iter, int *nit);

int forkal(int ip, int iq, int id, double *phi, double*theta, double *delta, int N, double *W, double *resid, int il, double *Y, double *AMSE);

double fcss_seas(double *b, int pq, void *params);

int css_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,
	double *phi, double *theta, double *PHI, double *THETA, double *wmean, double *var,double *loglik,double *hess);

double fas154_seas(double *b, int pq, void *params);

void checkroots(double *phi, int *p, double *theta, int *q, double *PHI, int *P, double *THETA, int *Q);

int as154_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,double *phi, double *theta, 
	double *PHI, double *THETA, double *wmean, double *var, double *loglik, double *hess, int cssml);

void free_alik_seas(alik_seas_object object);

void free_alik_css_seas(alik_css_seas_object object);

#ifdef __cplusplus
}
#endif

#endif /* FIT_H_ */
