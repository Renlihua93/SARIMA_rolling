// SPDX-License-Identifier: BSD-3-Clause
/* Internal declarations for seasonal ARIMA fitting (AS 154 / AS 182). */
#ifndef FIT_INTERNAL_H_
#define FIT_INTERNAL_H_

int as154_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,
	double *phi, double *theta, double *PHI, double *THETA,
	double *wmean, double *var, double *loglik, double *hess, int cssml);

int css_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,
	double *phi, double *theta, double *PHI, double *THETA,
	double *wmean, double *var, double *loglik, double *hess);

int forkal(int ip, int iq, int id, double *phi, double *theta, double *delta,
	int N, double *W, double *resid, int il, double *Y, double *AMSE);

void deld(int d, double *coef);
void delds(int D, int s, double *coef);
void conv(double *sig, int len_sig, double *filt, int len_filt, double *oup);

#endif /* FIT_INTERNAL_H_ */
