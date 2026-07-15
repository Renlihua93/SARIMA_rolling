// SPDX-License-Identifier: BSD-3-Clause
#ifndef SARIMA_H_
#define SARIMA_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle; internal layout is in sarima.c */
typedef struct SarimaState *sarima_object;

/* (p,d,q) x (P,D,Q)[s], series length N */
sarima_object sarima_init(int p, int d, int q, int s, int P, int D, int Q, int N);

void sarima_free(sarima_object obj);

/* 0 = MLE (default), 1 = CSS */
void sarima_setMethod(sarima_object obj, int method);

/* 0..7; default 5 (BFGS) */
void sarima_setOptMethod(sarima_object obj, int optmethod);

/* 1 = CSS starting values before MLE (default), 0 = skip */
void sarima_setCSSML(sarima_object obj, int cssml);

void sarima_exec(sarima_object obj, double *series);

void sarima_summary(sarima_object obj);

void sarima_vcov(sarima_object obj, double *vcov);

/* L-step forecast from full series inp; xpred and amse length L */
void sarima_predict(sarima_object obj, double *inp, int L, double *xpred, double *amse);

#ifdef __cplusplus
}
#endif

#endif /* SARIMA_H_ */
