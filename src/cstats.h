// SPDX-License-Identifier: BSD-3-Clause
/* cstats.h — 统计/分布/求根（SARIMA 精简） */
#ifndef CSTATS_H_
#define CSTATS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PIVAL 3.14159265358979323846264338327950288
#define XINFVAL 1.79e+308
#define XNINFVAL 2.2251e-308

#ifndef DBL_MAX_EXP
#define DBL_MAX_EXP +1024
#endif
#ifndef DBL_MIN_EXP
#define DBL_MIN_EXP -1021
#endif

double mean(double *vec, int N);

int cpoly(double *OPR, double *OPI, int DEGREE, double *ZEROR, double *ZEROI);
int polyroot(double *coeff, int DEGREE, double *ZEROR, double *ZEROI);

#ifdef __cplusplus
}
#endif

#endif /* CSTATS_H_ */
