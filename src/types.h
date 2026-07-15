// SPDX-License-Identifier: BSD-3-Clause
/*
 * types.h — 优化/NLS 回调类型（原 lnsrchmp.h 中的 typedef）
 * 供 linalg（NLS）与 optim 共用，避免头文件循环依赖。
 */
#ifndef TYPES_H_
#define TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct custom_function_set custom_function;

struct custom_function_set {
	double (*funcpt)(double *x, int N, void *params);
	void *params;
};

typedef struct custom_gradient_set custom_gradient;

struct custom_gradient_set {
	void (*funcgrad)(double *x, int N, double *g, void *params);
	void *params;
};

typedef struct custom_funcuni_set custom_funcuni;

struct custom_funcuni_set {
	double (*funcuni)(double x, void *params);
	void *params;
};

typedef struct custom_funcmult_set custom_funcmult;

struct custom_funcmult_set {
	void (*funcmult)(double *x, int M, int N, double *f, void *params);
	void *params;
};

typedef struct custom_jacobian_set custom_jacobian;

struct custom_jacobian_set {
	void (*jacobian)(double *x, int M, int N, double *jac, void *params);
	void *params;
};

#define FUNCPT_EVAL(F, x, N) (*((F)->funcpt))(x, N, (F)->params)
#define FUNCGRAD_EVAL(F, x, N, g) (*((F)->funcgrad))(x, N, (g), (F)->params)
#define FUNCUNI_EVAL(F, x) (*((F)->funcuni))(x, (F)->params)
#define FUNCMULT_EVAL(F, x, M, N, f) (*((F)->funcmult))(x, M, N, (f), (F)->params)
#define JACOBIAN_EVAL(F, x, M, N, jac) (*((F)->jacobian))(x, M, N, (jac), (F)->params)

#ifdef __cplusplus
}
#endif

#endif /* TYPES_H_ */
