// SPDX-License-Identifier: BSD-3-Clause
/*
 * linalg.c - SARIMA 矩阵运算
 *
 * 功能: 机器精度、向量范数、矩阵乘/转置、LU 分解与求逆。
 */

#include "linalg.h"


/*
 * macheps - 返回双精度机器 epsilon
 */
double macheps() {
	double macheps;
	macheps = 1.0;

	while ((macheps + 1.0) > 1.0) {
		macheps = macheps / 2.0;
	}

	macheps = macheps * 2;

	return macheps;
}


/*
 * signx - 实数符号 (-1 / +1)
 */
double signx(double x) {
	double sgn;
	if (x >= 0.) {
		sgn = 1.0;
	}
	else {
		sgn = -1.0;
	}
	return sgn;
}


/*
 * imax - 两整数取较大值
 */
int imax(int a, int b)
{
	return (a > b) ? a : b;
}


/*
 * imin - 两整数取较小值
 */
int imin(int a, int b)
{
	return (a < b) ? a : b;
}


/*
 * l2norm - 向量 L2 范数
 */
double l2norm(double *vec, int N) {
	double l2, sum;
	int i;
	sum = 0.;
	for (i = 0; i < N; ++i) {
		sum += vec[i] * vec[i];
	}
	l2 = sqrt(sum);
	return l2;
}


/*
 * array_max_abs - 向量元素绝对值的最大值
 */
double array_max_abs(double *array,int N) {
	int i;
	double m = 0.0;

	for (i = 0; i < N;++i) {
		if (fabs(array[i]) > m ) {
			m = fabs(array[i]);
		}
	}

	return m;
}


/*
 * scale - 矩阵标量乘法
 */
void scale(double *A, int rows, int cols, double alpha) {
	int N,i;
	 
	 
	N = rows * cols;
	
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (i = 0; i < N;++i) {
		A[i] = alpha * A[i];
	}
}


/*
 * madd - 矩阵相加
 */
void madd(double* A, double* B, double* C,int rows,int cols) {
	int N,i;
	 
	 
	N = rows * cols;
	
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (i = 0; i < N; ++i) {
		C[i] = A[i] + B[i];
	}
}


/*
 * msub - 矩阵相减
 */
void msub(double* A, double* B, double* C,int rows,int cols) {
	int N,i;
	 
	 
	N = rows * cols;
	
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (i = 0; i < N; ++i) {
		C[i] = A[i] - B[i];
	}
}


/*
 * nmult - 朴素矩阵乘法
 */
void nmult(double* A, double* B, double* C,int ra,int ca, int cb) {
	register int i,j,k;
	int u,v,t,rb;
	
	 
	 
	rb = ca;
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k,v,u,t)
#endif
	for (i = 0; i < ra; ++i) {
		for (j = 0; j < cb; ++j) {
			v = i * rb;
			u = i *cb;
			t = j + u;
			C[t] = 0.;
			for (k = 0; k < rb;++k) {
				C[t] += A[k + v] * B[j + k * cb];
			}
		}
	}


}


/*
 * mmult - 矩阵乘法 (调用 nmult)
 */
void mmult(double* A, double* B, double* C,int m,int n, int p) {
	(void)CUTOFF;
	nmult(A, B, C, m, n, p);
}


/*
 * stranspose - 矩阵转置 (标准实现)
 */
void stranspose(double *sig, int rows, int cols,double *col) {
	int t,u;
	register int i,j;
#ifdef _OPENMP
	#pragma omp parallel for private(i,j,t,u)
#endif
	for (i=0; i < rows; i++) {
		t = i * cols;
		u = 0;
		for (j=0; j < cols; j++) {
			col[u+i] = sig[j+t];
			u+=rows;
		}
	}
	
}


/*
 * mtranspose - 矩阵转置 (SARIMA 路径)
 */
void mtranspose(double *sig, int rows, int cols, double *col)
{
	stranspose(sig, rows, cols, col);
}


/*
 * pludecomp - 带部分主元的 LU 分解
 */
static int pludecomp(double *A,int N,int *ipiv) {
	int k,j,l,c1,c2,mind,tempi;
	double ld,mult,mval,temp;
	for(k=0;k < N;++k)
		ipiv[k] = k;
	
	for(k = 0; k < N-1; ++k) {
		
		mval = fabs(A[k*N + k]);
		mind = k;
		for (j=k+1; j < N;++j) {
			if (mval < fabs(A[j*N + k])) {
				mval = A[j*N + k];
				mind = j;
			}
		}
		
		if ( mind != k) {
			c1 = k *N;
			c2 = mind * N;
			tempi = ipiv[mind];
			ipiv[mind] = ipiv[k];
			ipiv[k] = tempi;
			for (j = 0; j < N;j++) {
				temp = A[c1 + j];
				*(A + c1 + j) = *(A + c2 + j);
				*(A + c2 + j) = temp;
			}
		}
		c2 = k*N;
		ld = A[c2 + k];
		if (ld != 0.) {
			for (j = k+1; j < N; ++j) {
				c1 = j*N;
				mult = A[c1+k] /= ld;
				
				for(l = k+1; l < N; ++l) {
					A[c1+l] -= mult * A[c2 + l];
				}
			}
		}
		
	}
	return 0;
	
}


/*
 * ludecomp - LU 分解入口
 */
void ludecomp(double *A,int N,int *ipiv) {
	pludecomp(A,N,ipiv);
}


/*
 * linsolve - LU 回代求解 Ax=b
 */
void linsolve(double *A,int N,double *b,int *ipiv,double *x) {
	int i,j,c1,l;
	double *y;
	double sum;
	
	y = (double*) malloc(sizeof(double) *N);
	 
	for(i = 0; i < N;++i) {
		y[i] = 0.;
		x[i] = 0.;
		if ( A[i*N + i] == 0.) {
			printf("Warning : The Matrix system does not have a unique solution");
			
		}
		
	}
	
	
	
	y[0] = b[ipiv[0]];
	for(i = 1; i < N; ++i) {
		sum = 0.;
		c1 = i*N;
		for(j = 0; j < i; ++j) {
			sum += y[j] * A[c1 + j];
		}
		y[i] = b[ipiv[i]] - sum;
	}
	
	
	
	x[N - 1] = y[N - 1]/A[N * N - 1];
	
	for (i = N - 2; i >= 0; i--) {
		sum = 0.;
		c1 = i*(N+1);
		l=0;
		for(j = i+1; j < N;j++) {
			l++;
			sum += A[c1 + l] * x[j];
		}
		x[i] = (y[i] - sum) / A[c1];
	}
	
	free(y);
}


/*
 * minverse - 矩阵求逆 (基于 LU)
 */
void minverse(double *A,int N,int *ipiv,double *inv) {
	int i,j,stride;
	double *col,*x;
	
	col = (double*) malloc(sizeof(double) * N);
	x = (double*) malloc(sizeof(double) * N);
	
	for (i = 0; i < N; ++i) {
		col[i] = 0.;
		x[i] = 0.;
	}
	
	for (i = 0; i < N; ++i) {
		col[i] = 1.;
		linsolve(A,N,col,ipiv,x);
		stride = i;
		for(j = 0; j < N;++j) {
			inv[stride] = x[j];
			stride+= N;
		}
		col[i] = 0.;
	}
		
	free(x);
	free(col);
}
