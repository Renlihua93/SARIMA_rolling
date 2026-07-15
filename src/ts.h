// SPDX-License-Identifier: BSD-3-Clause
/* ts.h — SARIMA 时序运算（精简） */
#ifndef TS_H_
#define TS_H_

#include "cstats.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/* (1-B)^d 与 (1-B^s)^D 算子系数 */
void deld(int d, double *C);
void delds(int D, int s, double *C);

/* 序列差分，返回差分后长度 */
int diff(double *sig, int N, int d, double *oup);
int diffs(double *sig, int N, int D, int s, double *oup);

void conv_direct(double *inp1, int N, double *inp2, int L, double *oup);
void conv(double *sig, int len_sig, double *filt, int len_filt, double *oup);

int archeck(int p, double *ar);
int invertroot(int q, double *ma);

#ifdef __cplusplus
}
#endif

#endif /* TS_H_ */
