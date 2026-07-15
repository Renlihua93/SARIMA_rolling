// SPDX-License-Identifier: BSD-3-Clause
/*
 * ts.c - SARIMA 时序工具
 *
 * 功能: 差分算子系数、序列差分、直接卷积、AR 平稳与 MA 可逆性检验。
 */

#include "ts.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


/*
 * poly_mult - 多项式系数相乘 (差分算子用)
 */
static int poly_mult(double *A, double *B, double *C, int lA, int lB)
{
	int lC, i, j, k;
	double temp;

	lC = lA + lB - 1;
	for (i = 0; i < lC; ++i)
		C[i] = 0.0;
	for (i = 0; i < lC; ++i) {
		temp = 0.0;
		for (j = 0; j < lA; ++j) {
			for (k = 0; k < lB; ++k) {
				if (j + k == i)
					temp += A[j] * B[k];
			}
		}
		C[i] = temp;
	}
	return lC;
}


/*
 * deld - 构造普通差分算子系数
 */
void deld(int d, double *C)
{
	int i, j;
	double *vec, *oup;

	if (d <= 0) {
		C[0] = 1.0;
		return;
	}
	vec = (double *)malloc(2 * sizeof(double));
	oup = (double *)malloc((size_t)(d + 1) * sizeof(double));
	vec[0] = 1.0;
	vec[1] = -1.0;
	oup[0] = 1.0;
	for (i = 0; i < d; ++i) {
		poly_mult(oup, vec, C, i + 1, 2);
		for (j = 0; j < i + 2; ++j)
			oup[j] = C[j];
	}
	free(vec);
	free(oup);
}


/*
 * delds - 构造季节差分算子系数
 */
void delds(int D, int s, double *C)
{
	int i, j;
	double *vec, *oup;

	if (D <= 0) {
		C[0] = 1.0;
		return;
	}
	vec = (double *)calloc((size_t)(s + 1), sizeof(double));
	oup = (double *)malloc((size_t)(D * s + 1) * sizeof(double));
	vec[0] = 1.0;
	vec[s] = -1.0;
	oup[0] = 1.0;
	for (i = 0; i < D; ++i) {
		poly_mult(oup, vec, C, i * s + 1, s + 1);
		for (j = 0; j < s * (i + 1) + 1; ++j)
			oup[j] = C[j];
	}
	free(vec);
	free(oup);
}


/*
 * diff - 对序列施加普通差分 (d 阶)
 */
int diff(double *sig, int N, int d, double *oup)
{
	int Noup, i, j;
	double *coeff, sum;

	if (d <= 0) {
		for (i = 0; i < N; ++i)
			oup[i] = sig[i];
		return N;
	}
	coeff = (double *)malloc((size_t)(d + 1) * sizeof(double));
	deld(d, coeff);
	Noup = N - d;
	for (i = d; i < N; ++i) {
		sum = 0.0;
		for (j = 1; j < d + 1; ++j)
			sum += sig[i - j] * coeff[j];
		oup[i - d] = sum + coeff[0] * sig[i];
	}
	free(coeff);
	return Noup;
}


/*
 * diffs - 对序列施加季节差分 (D 阶, 周期 s)
 */
int diffs(double *sig, int N, int D, int s, double *oup)
{
	int Noup, i, j, d;
	double *coeff, sum;

	if (D <= 0) {
		for (i = 0; i < N; ++i)
			oup[i] = sig[i];
		return N;
	}
	d = D * s;
	coeff = (double *)malloc((size_t)(d + 1) * sizeof(double));
	delds(D, s, coeff);
	Noup = N - d;
	for (i = d; i < N; ++i) {
		sum = 0.0;
		for (j = s; j < d + 1; j += s)
			sum += sig[i - j] * coeff[j];
		oup[i - d] = sum + coeff[0] * sig[i];
	}
	free(coeff);
	return Noup;
}


/*
 * conv_direct - 直接法一维卷积
 */
void conv_direct(double *inp1, int N, double *inp2, int L, double *oup)
{
	int M, k, m, i;
	double t1, tmin;

	M = N + L - 1;
	i = 0;
	if (N >= L) {
		for (k = 0; k < L; k++) {
			oup[k] = 0.0;
			for (m = 0; m <= k; m++)
				oup[k] += inp1[m] * inp2[k - m];
		}
		for (k = L; k < M; k++) {
			oup[k] = 0.0;
			i++;
			t1 = (double)L + i;
			tmin = MIN(t1, N);
			for (m = i; m < tmin; m++)
				oup[k] += inp1[m] * inp2[k - m];
		}
	} else {
		for (k = 0; k < N; k++) {
			oup[k] = 0.0;
			for (m = 0; m <= k; m++)
				oup[k] += inp2[m] * inp1[k - m];
		}
		for (k = N; k < M; k++) {
			oup[k] = 0.0;
			i++;
			t1 = (double)N + i;
			tmin = MIN(t1, L);
			for (m = i; m < tmin; m++)
				oup[k] += inp2[m] * inp1[k - m];
		}
	}
}


/*
 * conv - 一维实序列卷积
 */
void conv(double *sig, int len_sig, double *filt, int len_filt, double *oup)
{
	if (len_sig >= len_filt)
		conv_direct(sig, len_sig, filt, len_filt, oup);
	else
		conv_direct(filt, len_filt, sig, len_sig, oup);
}


/*
 * archeck - 检验 AR 特征多项式根是否在单位圆内 (平稳性)
 */
int archeck(int p, double *ar)
{
	int check, N, i;
	double *coeff, *zeror, *zeroi, mod;

	if (p <= 0)
		return 1;
	N = p + 1;
	check = 1;
	coeff = (double *)malloc((size_t)N * sizeof(double));
	zeror = (double *)malloc((size_t)p * sizeof(double));
	zeroi = (double *)malloc((size_t)p * sizeof(double));
	coeff[0] = 1.0;
	for (i = 0; i < p; ++i)
		coeff[i + 1] = -ar[i];
	polyroot(coeff, p, zeror, zeroi);
	for (i = 0; i < p; ++i) {
		mod = sqrt(zeror[i] * zeror[i] + zeroi[i] * zeroi[i]);
		if (mod < 1.0) {
			check = 0;
			break;
		}
	}
	free(coeff);
	free(zeror);
	free(zeroi);
	return check;
}


/*
 * invertroot - 检验 MA 特征多项式根是否在单位圆外 (可逆性)
 */
int invertroot(int q, double *ma)
{
	int i, index, retval, fail, rcheck, qn, i1, j;
	double *temp, *zeror, *zeroi, *xr, *xi, *yr, *yi;
	int *ind;
	double mod, tempr, tempi;

	retval = 0;
	index = -1;
	for (i = 0; i < q; ++i) {
		if (ma[i] != 0.0)
			index = i;
	}
	if (index == -1)
		return retval;

	temp = (double *)malloc((size_t)(q + 1) * sizeof(double));
	zeror = (double *)malloc((size_t)q * sizeof(double));
	zeroi = (double *)malloc((size_t)q * sizeof(double));
	ind = (int *)malloc((size_t)q * sizeof(int));
	index++;
	temp[0] = 1.0;
	for (i = 1; i <= index; ++i)
		temp[i] = ma[i - 1];
	qn = index;
	xr = (double *)malloc((size_t)(qn + 1) * sizeof(double));
	xi = (double *)malloc((size_t)(qn + 1) * sizeof(double));
	yr = (double *)malloc((size_t)(qn + 1) * sizeof(double));
	yi = (double *)malloc((size_t)(qn + 1) * sizeof(double));
	fail = polyroot(temp, qn, zeror, zeroi);
	if (fail == 1) {
		free(zeror);
		free(zeroi);
		free(temp);
		free(ind);
		free(xr);
		free(xi);
		free(yr);
		free(yi);
		return retval;
	}
	rcheck = 0;
	for (i = 0; i < qn; ++i) {
		mod = zeror[i] * zeror[i] + zeroi[i] * zeroi[i];
		ind[i] = (mod < 1.0) ? 1 : 0;
		if (ind[i])
			rcheck++;
	}
	if (rcheck == 0) {
		free(zeror);
		free(zeroi);
		free(temp);
		free(ind);
		free(xr);
		free(xi);
		free(yr);
		free(yi);
		return retval;
	}
	retval = 1;
	if (index == 1) {
		ma[0] = 1.0 / ma[0];
		for (i = 1; i < q; ++i)
			ma[i] = 0.0;
		free(zeror);
		free(zeroi);
		free(temp);
		free(ind);
		free(xr);
		free(xi);
		free(yr);
		free(yi);
		return retval;
	}
	for (i = 0; i < index; ++i) {
		if (ind[i] == 1) {
			mod = zeror[i] * zeror[i] + zeroi[i] * zeroi[i];
			zeror[i] /= mod;
			zeroi[i] = -zeroi[i] / mod;
		}
	}
	xr[0] = 1.0;
	xi[0] = 0.0;
	for (i = 0; i < qn; ++i) {
		i1 = i + 1;
		mod = zeror[i] * zeror[i] + zeroi[i] * zeroi[i];
		tempr = zeror[i] / mod;
		tempi = -zeroi[i] / mod;
		xr[i1] = xi[i1] = 0.0;
		for (j = 1; j <= i1; ++j) {
			yr[j] = tempr * xr[j - 1] - tempi * xi[j - 1];
			yi[j] = tempr * xi[j - 1] + tempi * xr[j - 1];
			yr[j] = xr[j] - yr[j];
			yi[j] = xi[j] - yi[j];
		}
		for (j = 1; j <= i1; ++j) {
			xr[j] = yr[j];
			xi[j] = yi[j];
		}
	}
	for (i = 0; i < qn; ++i)
		ma[i] = xr[i + 1];
	for (i = qn; i < q; ++i)
		ma[i] = 0.0;
	free(zeror);
	free(zeroi);
	free(temp);
	free(ind);
	free(xr);
	free(xi);
	free(yr);
	free(yi);
	return retval;
}
