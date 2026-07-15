// SPDX-License-Identifier: BSD-3-Clause
/*
 * sarima.c - SARIMA 对外 API
 *
 * 功能: 模型创建、配置、拟合、预测、摘要打印、协方差导出与释放。
 * 依赖 fit(似然) 与 ts(差分、Kalman 预测)。
 */

#include "sarima.h"
#include "fit_internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



struct SarimaState {
	int N, Nused;
	int method, optmethod, cssml;
	int p, d, q, s, P, D, Q;
	int M, ncoeff, lvcov, retval;
	double mean, var, loglik, aic;
	double *phi, *theta, *PHI, *THETA, *res, *vcov;
	double params[];
};




/*
 * opt_method_name - 将优化器编号映射为可读名称
 */
static const char *opt_method_name(int m)
{
	static const char *names[] = {
		"Nelder-Mead", "Newton Line Search",
		"Newton Trust Region - Hook Step", "Newton Trust Region - Double Dog-Leg",
		"Conjugate Gradient", "BFGS", "L-BFGS", "BFGS More-Thuente Line Search"
	};
	if (m >= 0 && m <= 7)
		return names[m];
	return "Unknown";
}




/*
 * exit_message - 将优化返回码映射为可读说明
 */
static const char *exit_message(int code)
{
	switch (code) {
	case 0: return "Input Error";
	case 1: return "Probable Success";
	case 4: return "Optimization Routine didn't converge";
	case 15: return "Optimization Routine Encountered Inf/Nan Values";
	default: return "Unknown";
	}
}






/*
 * sarima_init - 分配并初始化 SARIMA 模型对象
 */
sarima_object sarima_init(int p, int d, int q, int s, int P, int D, int Q, int N)
{
	sarima_object obj;
	int i, npar, nres, M;

	if (p < 0 || d < 0 || q < 0 || N <= 0 || P < 0 || D < 0 || Q < 0 || s < 0) {
		fprintf(stderr, "sarima_init: orders and N must be non-negative; N > 0.\n");
		exit(-1);
	}

	M = (d == 0 && D == 0) ? 1 : 0;
	npar = p + q + P + Q;
	nres = N - d - s * D;

	obj = (sarima_object)malloc(sizeof(struct SarimaState)
		+ (size_t)npar * sizeof(double)
		+ (size_t)nres * sizeof(double)
		+ (size_t)(npar + M) * (npar + M) * sizeof(double));
	if (!obj)
		return NULL;

	memset(obj, 0, sizeof(struct SarimaState));
	obj->p = p; obj->d = d; obj->q = q;
	obj->s = s; obj->P = P; obj->D = D; obj->Q = Q;
	obj->N = N; obj->Nused = nres;
	obj->M = M; obj->ncoeff = npar + M;
	obj->lvcov = obj->ncoeff * obj->ncoeff;
	obj->method = 0;
	obj->optmethod = 5;
	obj->cssml = 1;
	obj->var = 1.0;

	obj->phi = obj->params;
	obj->theta = obj->params + p;
	obj->PHI = obj->params + p + q;
	obj->THETA = obj->params + p + q + P;
	obj->res = obj->params + npar;
	obj->vcov = obj->res + nres;

	for (i = 0; i < npar; ++i)
		obj->params[i] = 0.0;

	return obj;
}






/*
 * sarima_setMethod - 设置估计方法: 0 为 MLE, 1 为 CSS
 */
void sarima_setMethod(sarima_object obj, int value)
{
	if (value == 0 || value == 1)
		obj->method = value;
	else
		fprintf(stderr, "sarima_setMethod: use 0 (MLE) or 1 (CSS).\n");
}






/*
 * sarima_setCSSML - 设置 MLE 是否从 CSS 初值启动
 */
void sarima_setCSSML(sarima_object obj, int cssml)
{
	if (cssml == 0 || cssml == 1)
		obj->cssml = cssml;
	else {
		fprintf(stderr, "sarima_setCSSML: use 0 or 1.\n");
		exit(-1);
	}
}






/*
 * sarima_setOptMethod - 设置优化器编号 0..7
 */
void sarima_setOptMethod(sarima_object obj, int value)
{
	if (value >= 0 && value <= 7)
		obj->optmethod = value;
	else
		fprintf(stderr, "sarima_setOptMethod: use 0..7.\n");
}






/*
 * sarima_exec - 对输入序列执行季节 ARIMA 拟合
 */
void sarima_exec(sarima_object obj, double *inp)
{
	if (obj->method == 0) {
		obj->retval = as154_seas(inp, obj->N, obj->optmethod,
			obj->p, obj->d, obj->q, obj->s, obj->P, obj->D, obj->Q,
			obj->phi, obj->theta, obj->PHI, obj->THETA,
			&obj->mean, &obj->var, &obj->loglik, obj->vcov, obj->cssml);
		obj->loglik = -0.5 * (obj->Nused * (2.0 * obj->loglik + 1.0 + log(2.0 * M_PI)));
		obj->aic = -2.0 * obj->loglik + 2.0 * obj->ncoeff + 2.0;
	} else if (obj->method == 1) {
		obj->retval = css_seas(inp, obj->N, obj->optmethod,
			obj->p, obj->d, obj->q, obj->s, obj->P, obj->D, obj->Q,
			obj->phi, obj->theta, obj->PHI, obj->THETA,
			&obj->mean, &obj->var, &obj->loglik, obj->vcov);
		obj->loglik = -0.5 * (obj->Nused * (2.0 * obj->loglik + 1.0 + log(2.0 * M_PI)));
	} else {
		fprintf(stderr, "sarima_exec: unknown method %d.\n", obj->method);
	}
}







/*
 * expand_sarima_operators - 季节 AR/MA 展开为非季节算子
 */
static void expand_sarima_operators(sarima_object obj, int ir,
	double *phi, double *theta)
{
	int i, j, p, q, P, Q, s;

	p = obj->p; q = obj->q; P = obj->P; Q = obj->Q; s = obj->s;
	for (i = 0; i < ir; ++i)
		phi[i] = theta[i] = 0.0;

	for (i = 0; i < p; ++i)
		phi[i] = obj->phi[i];
	for (i = 0; i < q; ++i)
		theta[i] = -obj->theta[i];

	for (j = 0; j < P; ++j) {
		phi[(j + 1) * s - 1] += obj->PHI[j];
		for (i = 0; i < p; ++i)
			phi[(j + 1) * s + i] -= obj->phi[i] * obj->PHI[j];
	}
	for (j = 0; j < Q; ++j) {
		theta[(j + 1) * s - 1] -= obj->THETA[j];
		for (i = 0; i < q; ++i)
			theta[(j + 1) * s + i] += obj->theta[i] * obj->THETA[j];
	}
}






/*
 * sarima_predict - L 步超前预测并给出均方误差
 */
void sarima_predict(sarima_object obj, double *inp, int L, double *xpred, double *amse)
{
	int N, d, D, s, p, q, P, Q, ip, iq, ir, id, i;
	double *coef1, *coef2, *delta, *W, *resid, *phi, *theta;
	double wmean;

	N = obj->N; d = obj->d; D = obj->D; s = obj->s;
	p = obj->p; q = obj->q; P = obj->P; Q = obj->Q;

	ip = p + s * P;
	iq = q + s * Q;
	ir = ip;
	if (ir < 1 + q + s * Q)
		ir = 1 + q + s * Q;
	id = d + D * s;

	coef1 = (double *)malloc((size_t)(d + 1) * sizeof(double));
	coef2 = (double *)malloc((size_t)(D * s + 1) * sizeof(double));
	delta = (double *)malloc((size_t)(d + D * s + 1) * sizeof(double));
	W = (double *)malloc((size_t)N * sizeof(double));
	resid = (double *)malloc((size_t)N * sizeof(double));
	phi = (double *)malloc((size_t)ir * sizeof(double));
	theta = (double *)malloc((size_t)ir * sizeof(double));

	coef1[0] = coef2[0] = 1.0;
	wmean = (d == 0 && D == 0) ? obj->mean : 0.0;

	if (d > 0)
		deld(d, coef1);
	if (D > 0)
		delds(D, s, coef2);
	conv(coef1, d + 1, coef2, D * s + 1, delta);

	for (i = 1; i <= d + D * s; ++i)
		delta[i] = -delta[i];

	for (i = 0; i < N; ++i) {
		W[i] = inp[i] - wmean;
		resid[i] = obj->res[i];
	}

	expand_sarima_operators(obj, ir, phi, theta);
	forkal(ip, iq, id, phi, theta, delta + 1, N, W, resid, L, xpred, amse);

	for (i = 0; i < L; ++i)
		xpred[i] += wmean;

	free(coef1); free(coef2); free(delta);
	free(W); free(resid); free(phi); free(theta);
}






/*
 * sarima_vcov - 复制参数协方差矩阵
 */
void sarima_vcov(sarima_object obj, double *vcov)
{
	memcpy(vcov, obj->vcov, (size_t)obj->lvcov * sizeof(double));
}






/*
 * sarima_summary - 打印系数、标准误、对数似然与 AIC
 */
void sarima_summary(sarima_object obj)
{
	int i, pq, t;

	pq = obj->ncoeff;
	printf("\n\nExit Status \n");
	printf("Return Code : %d \n", obj->retval);
	printf("Exit Message : %s\n", exit_message(obj->retval));

	printf("\n\nARIMA Seasonal Order : (%d, %d, %d) * (%d, %d, %d)[%d]\n\n",
		obj->p, obj->d, obj->q, obj->P, obj->D, obj->Q, obj->s);

	printf("%-20s%-20s%-20s \n\n", "Coefficients", "Value", "Standard Error");
	for (i = 0; i < obj->p; ++i)
		printf("AR%-15d%-20g%-20g \n", i + 1, obj->phi[i], sqrt(obj->vcov[i + pq * i]));
	for (i = 0; i < obj->q; ++i) {
		t = obj->p + i;
		printf("MA%-15d%-20g%-20g \n", i + 1, obj->theta[i], sqrt(obj->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->P; ++i) {
		t = obj->p + obj->q + i;
		printf("SAR%-14d%-20g%-20g \n", i + 1, obj->PHI[i], sqrt(obj->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->Q; ++i) {
		t = obj->p + obj->q + obj->P + i;
		printf("SMA%-14d%-20g%-20g \n", i + 1, obj->THETA[i], sqrt(obj->vcov[t + pq * t]));
	}
	t = obj->p + obj->q + obj->P + obj->Q;
	if (obj->M == 1)
		printf("%-17s%-20g%-20g \n", "MEAN", obj->mean, sqrt(obj->vcov[t + pq * t]));
	else
		printf("%-17s%-20g \n", "MEAN", obj->mean);

	printf("\n%-17s%-20g \n", "SIGMA^2", obj->var);
	printf("\nESTIMATION METHOD : %s\n\n",
		obj->method == 0 ? "MLE" : "CSS");
	printf("OPTIMIZATION METHOD : %s\n\n", opt_method_name(obj->optmethod));

	if (obj->method == 0 || obj->method == 1)
		printf("Log Likelihood : %g \n\n", obj->loglik);
	else
		printf("Log Likelihood : Unavailable \n\n");

	if (obj->method == 0)
		printf("AIC criterion : %g \n\n", obj->aic);
	else
		printf("AIC Criterion : Unavailable \n\n");
}






/*
 * sarima_free - 释放 SARIMA 模型对象
 */
void sarima_free(sarima_object obj)
{
	free(obj);
}
