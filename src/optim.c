// SPDX-License-Identifier: BSD-3-Clause
/*
 * optim.c - SARIMA 数值优化
 *
 * 功能: 有限差分梯度、BFGS 无约束极小化、fminunc 入口、Hessian 数值近似。
 */

#include "optim.h"


/*
 * grad_fd - 梯度: 有限差分或用户提供的梯度回调
 */
int grad_fd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
		double eps2, double *f) {
	int retval;
	retval = 0;
	if (funcgrad == NULL) {
		
		retval = grad_calc(funcpt,x,N,dx,eps2,f);
	} else {
		
		FUNCGRAD_EVAL(funcgrad,x,N,f);
	}

	return retval;

}


/*
 * grad_calc - 数值梯度 (中心差分)
 */
int grad_calc(custom_function *funcpt, double *x, int N, double *dx, double eps2, double *f) {
	int i, j,retval;
	double step, fd, stepmax;
	double *xi;

	fd = eps2; 
	retval = 0;
	xi = (double*)malloc(sizeof(double)*N);

	for (i = 0; i < N; ++i) {
		if (fabs(x[i]) >= 1.0 / fabs(dx[i])) {
			stepmax = x[i];
		}
		else {
			stepmax = signx(x[i]) * 1.0 / fabs(dx[i]);
		}
		step = fd * stepmax;
		for (j = 0; j < N; ++j) {
			xi[j] = x[j];
		}
		xi[i] += step;
		f[i] = (FUNCPT_EVAL(funcpt, xi, N) - FUNCPT_EVAL(funcpt, x, N)) / step;
		if (f[i] >= DBL_MAX || f[i] <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(xi);
			return 15;
		}
		if (f[i] != f[i]) {
			printf("Program Exiting as the function returns NaN");
			free(xi);
			return 15;
		}
		
	}

	free(xi);
	return retval;

}


/*
 * linsolve_lower - 下三角 (Cholesky 因子) 线性方程求解
 */
void linsolve_lower(double *L,int N,double *b,double *x) {
	int i,j,c1,l;
	double *y,*A;
	double sum;
	
	y = (double*) malloc(sizeof(double) *N);
	A = (double*) malloc(sizeof(double) *N *N);
	
	for(i = 0; i < N;++i) {
		y[i] = 0.;
		x[i] = 0.;
		if ( L[i*N + i] == 0.) {
			printf("The Matrix system does not have a unique solution");
			exit(1);
		}
		
	}
	
	
	
	y[0] = b[0]/L[0];
	for(i = 1; i < N; ++i) {
		sum = 0.;
		c1 = i*N;
		for(j = 0; j < i; ++j) {
			sum += y[j] * L[c1 + j];
		}
		y[i] = (b[i] - sum)/L[c1+i];
	}
	
	mtranspose(L,N,N,A);
	
	
	
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
	free(A);
}


/*
 * hessian_fd - Hessian 有限差分近似
 */
int hessian_fd(custom_function *funcpt,double *x,int N,double *dx,double eps,double *f) {
	int i,j,retval;
	double stepi,stepj,fd,stepmax;
	double *xi,*xj,*xij;
	
	retval = 0;

	fd = pow((double) eps,1.0/3.0);
	xi = (double*) malloc(sizeof(double) *N);
	xj = (double*) malloc(sizeof(double) *N);
	xij = (double*) malloc(sizeof(double) *N);
	
	for(i = 0; i < N;++i) {
		if (fabs(x[i]) >= 1.0 / fabs(dx[i])) {
			stepmax = x[i];
		} else {
			stepmax = signx(x[i]) * 1.0 / fabs(dx[i]);
		}
		stepi = fd * stepmax;
		for(j = 0; j < N;++j) {
			xi[j] = x[j];
			xj[j]= x[j];
			xij[j] = x[j];
		}
		xi[i] += stepi;
		xij[i] += stepi;
		for(j = 0; j < N;++j) {
			if (fabs(x[j]) >= 1.0 / fabs(dx[j])) {
				stepmax = x[j];
			} else {
				stepmax = signx(x[j]) * 1.0 / fabs(dx[j]);
			}
			stepj = fd * stepmax;
			xj[j] += stepj;
			xij[j] += stepj;

			f[i*N+j] = ((FUNCPT_EVAL(funcpt,xij,N) - FUNCPT_EVAL(funcpt,xi,N)) - (FUNCPT_EVAL(funcpt,xj,N) - FUNCPT_EVAL(funcpt,x,N)))/(stepi * stepj);
			if (f[i*N+j] >= DBL_MAX || f[i*N+j] <= -DBL_MAX) {
				printf("Program Exiting as the function value exceeds the maximum double value");
				free(xi);
				free(xj);
				free(xij);
				return 15;
			}
			if (f[i*N+j] != f[i*N+j]) {
				printf("Program Exiting as the function returns NaN");
				free(xi);
				free(xj);
				free(xij);
				return 15;
			}
			xj[j] -= stepj;
			xij[j] -= stepj;
		}

	}
	
	free(xi);
	free(xj);
	free(xij);
	return retval;
	
}


/*
 * lnsrch - 线搜索求可接受步长
 */
int lnsrch(custom_function *funcpt,double *xi,double *jac,double *p,int N,double * dx,double maxstep,double stol,double *x) {
	int retval,i,iter,MAXITER;
	double alpha,lambda,lambdamin,funcf,funci,lambdaprev,lambdatemp,funcprev;
	double lambda2,lambdaprev2,ll,den,rell,nlen;
	double *slopei,*temp1,*temp2,*ab,*rcheck,*pl;
	
	slopei = (double*) malloc(sizeof(double) *1);
	temp1 = (double*) malloc(sizeof(double) *4);
	temp2 = (double*) malloc(sizeof(double) *2);
	ab = (double*) malloc(sizeof(double) *2);
	rcheck = (double*) malloc(sizeof(double) *N);
	pl = (double*) malloc(sizeof(double) *N);
	retval = 100;
	alpha = 1e-04;
	lambda = 1.0;
	nlen = 0.0;
	MAXITER = 1000;
	funcprev = 1.0; 
	lambdaprev = 1.0; 
	for(i = 0; i < N;++i) {
		nlen += dx[i] * p[i] * dx[i] * p[i];
	}
	nlen = sqrt(nlen);
	iter = 0;
	if (nlen > maxstep) {
		scale(p,1,N,maxstep/nlen);
		nlen = maxstep;
	}
	
	mmult(jac,p,slopei,1,N,1);
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 /fabs(dx[i])) {
			den = fabs(xi[i]);
		} else {
			den = 1.0 /fabs(dx[i]);
		}
		rcheck[i] = p[i]/den;
	}
	
	rell = array_max_abs(rcheck,N);
	
	lambdamin = stol/rell;
	
	
	funci = FUNCPT_EVAL(funcpt, xi, N);
	if (funci >= DBL_MAX || funci <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		free(slopei);
		free(temp1);
		free(temp2);
		free(ab);
		free(rcheck);
		free(pl);
		return 15;
	}
	if (funci != funci) {
		printf("Program Exiting as the function returns NaN");
		free(slopei);
		free(temp1);
		free(temp2);
		free(ab);
		free(rcheck);
		free(pl);
		return 15;
	}
	while (retval > 1 && iter < MAXITER) {
		iter++;
		for(i = 0; i < N;++i) {
			pl[i] = p[i] * lambda;
		}
		madd(xi,pl,x,1,N);
		funcf = FUNCPT_EVAL(funcpt, x, N);
		if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(slopei);
			free(temp1);
			free(temp2);
			free(ab);
			free(rcheck);
			free(pl);
			return 15;
		}
		if (funcf != funcf) {
			printf("Program Exiting as the function returns NaN");
			free(slopei);
			free(temp1);
			free(temp2);
			free(ab);
			free(rcheck);
			free(pl);
			return 15;
		}

		
		
		if (funcf <= funci + alpha *lambda *slopei[0]) {
			retval = 0;
		} else if (lambda < lambdamin) {
			retval = 1;
			for (i = 0; i < N;++i) {
				x[i] = xi[i]; 
			}
		} else {
			if (lambda == 1.0) {
				lambdatemp = - slopei[0] / (2.0 * (funcf - funci - slopei[0])); 
			} else {
				lambda2 = lambda * lambda;
				lambdaprev2 = lambdaprev * lambdaprev;
				ll = lambda - lambdaprev;
				temp1[0] = 1.0 / lambda2; temp1[1] = -1.0 /lambdaprev2;
				temp1[2] = - lambdaprev / lambda2; temp1[3] = lambda /lambdaprev2;
				temp2[0] = funcf - funci - lambda * slopei[0];
				temp2[1] = funcprev - funci - lambdaprev * slopei[0];
				mmult(temp1,temp2,ab,2,2,1);
				scale(ab,1,2,1.0/ll);
				if (ab[0] == 0.0) {
					lambdatemp = - slopei[0] / (2.0 * ab[1]);
				} else {
					lambdatemp = (-ab[1] + sqrt( ab[1] * ab[1] - 3.0 * ab[0] *slopei[0]))/ (3.0 * ab[0]);
				}
				
				if (lambdatemp > 0.5 * lambda) {
					lambdatemp = 0.5 * lambda;
				}
			}
			lambdaprev = lambda;
			funcprev = funcf;
			if (lambdatemp <= 0.1 * lambda) {
				lambda = 0.1 * lambda;
			} else {
				lambda = lambdatemp;
			}
		}
	
	}
	
	free(slopei);
	free(temp1);
	free(temp2);
	free(ab);
	free(rcheck);
	free(pl);
	return retval;
}


/*
 * stopcheck - 优化停止判断 (梯度/步长/函数值)
 */
int stopcheck(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol,int retval) {
	int rcode,i;
	double num,den;
	double stop0;
	double *scheck;
	
	rcode = 0;	
	if (retval == 1) {
		rcode = 3;
		return rcode;
	}
	if (retval == 5) {
		rcode = 5;
		return rcode;
	}
	if (retval == 15) {
		rcode = 15;
		return rcode;
	}
	scheck = (double*) malloc(sizeof(double) *N);
	
	if (fabs(fx) > fabs(fsval)) {
			den = fabs(fx);
	} else {
			den = fabs(fsval);
	}
	for(i = 0; i < N;++i) {
		if (fabs(xf[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xf[i]);
		} else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}
	
	stop0 = array_max_abs(scheck,N);

	if (stop0 <= gtol) {
		rcode = 1;
	} else {
		for(i = 0; i < N;++i) {
			if (fabs(xf[i]) > 1.0 / fabs(dx[i])) {
				den = fabs(xf[i]);
			} else {
				den = 1.0 / fabs(dx[i]);
			}
			num = fabs(xf[i] - xc[i]);
			scheck[i] = num / den;
		}
		stop0 = array_max_abs(scheck,N);
		if (stop0 <= stol) {
			rcode = 2;
		}
	}
	
	free(scheck);
	return rcode;
}


/*
 * jrotate - Givens 旋转 (QR 秩1更新辅助)
 */
static void jrotate(double *A,int N,double a,double b,int i) {
	int j,r,r1;
	double c,s,den,y,w;
	r = i * N;
	r1 = r + N;
	if (a == 0.0) {
		c = 0.0;
		s = signx(b);
	} else {
		den = sqrt(a*a+b*b);
		c = a / den;
		s = b / den;
	}
	
	for (j = i; j < N;++j) {
		y = A[r+j];
		w = A[r1+j];
		A[r+j] = c*y - s*w;
		A[r1+j] = s*y + c*w;
	}
}


/*
 * qrupdate - QR 因子秩1更新
 */
static void qrupdate(double *R,int N,double *u,double *v) {
	int i,j,k,l;
	
	for(i = 1; i < N;++i) {
		R[i*N+i-1] = 0;
	}
	
	k = N-1;
	while (u[k] == 0 && k > 0) {
		k--;
	}
	
	for (i = k-1; i >= 0;--i) {
		jrotate(R,N,u[i],-u[i+1],i);
		if (u[i] == 0) {
			u[i] = fabs(u[i+1]);
		} else {
			u[i] = sqrt(u[i]*u[i] + u[i+1]*u[i+1]);
		}
	}
	
	for(i = 0; i < N;++i) {
		R[i] += u[0] * v[i];
	}
	
	for(i = 0; i < k;++i) {
		j = i * N;
		l = j + N;
		jrotate(R,N,R[j+i],-R[l+i],i);
	}
}


/*
 * inithess_lower - BFGS 初始逆 Hessian (下三角因子)
 */
static void inithess_lower(double *L,int N,double fi,double fsval,double *dx) {
	int i,j,ct;
	double temp;
	
	if (fabs(fi) > fsval) {
		temp = fabs(fi);
	} else {
		temp = fsval;
	}
	
	temp = sqrt(temp);
	
	for(i = 0; i < N;++i) {
		ct = i *N;
		L[ct+i] = temp * dx[i];
		for(j = 0; j < i;++j) {
			L[ct+j] = 0;
		}
	}
	
}


/*
 * bfgs_factored - BFGS 逆 Hessian 因子更新
 */
void bfgs_factored(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf) {
	int i,j,supd,ct;
	double sn,yn,fd,yt,jacm,alpha,temp3;
	double *sk,*yk,*temp,*temp2,*t,*u;
	
	sk = (double*) malloc(sizeof(double) *N);
	yk = (double*) malloc(sizeof(double) *N);
	temp = (double*) malloc(sizeof(double) *1);
	temp2 = (double*) malloc(sizeof(double) *1);
	t = (double*) malloc(sizeof(double) *N);
	u = (double*) malloc(sizeof(double) *N);
	
	msub(xf,xi,sk,1,N);
	msub(jacf,jac,yk,1,N);
	
	sn = l2norm(sk,N);
	yn = l2norm(yk,N);
	fd = sqrt(eps);
	
	mmult(yk,sk,temp,1,N,1);
	
	if (temp[0] >= fd*sn*yn) {
		for(i = 0; i < N;++i) {
			t[i] = 0.0;
			for(j = i; j < N;++j) {
				ct = j * N;
				t[i] += H[ct + i] * sk[j];
			}
		}
		mmult(t,t,temp2,1,N,1);
		alpha = sqrt(temp[0]/temp2[0]);
		supd = 1;
		for(i = 0; i < N;++i) {
			temp3 = 0.0;
			ct = i * N;
			for(j = 0; j < i+1;++j) {
				temp3 += H[ct + j] * t[j];
			}
			yt = fabs(yk[i] - temp3);
			if (jac[i] > jacf[i]) {
				jacm = jac[i];
			} else {
				jacm = jacf[i];
			}
			if (yt >= fd * jacm) {
				supd = 0;
			}
			u[i] = yk[i] - alpha*temp3;
		}
		if (supd == 0) {
			temp3 = 1.0 / sqrt(temp[0] * temp2[0]);
			for(i = 0; i < N; ++i) {
				t[i] *= temp3;
			}
			for(i = 1; i < N; ++i) {
				ct = i *N;
				for(j = 0; j < i;++j) {
					H[j*N+i] = H[ct+j];
				}
			}
			qrupdate(H,N,t,u);
			for(i = 1; i < N; ++i) {
				ct = i *N;
				for(j = 0; j < i;++j) {
					H[ct+j] = H[j*N+i];
				}
			}
		} 
	} 
	
	free(sk);
	free(yk);
	free(temp);
	free(temp2);
	free(t);
	free(u);
}


/*
 * bfgs_min - BFGS 无约束极小化主循环
 */
int bfgs_min(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval,double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double stol,double *xf)  {
	int rcode,gfdcode;
	int i,siter,retval;
	double dt1,dt2;
	double fx,num,den,stop0,fxf,eps2;
	double *jac,*hess,*scheck,*xc,*L,*step,*jacf;
	
	jac = (double*) malloc(sizeof(double) *N);
	scheck = (double*) malloc(sizeof(double) *N);
	xc = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	jacf = (double*) malloc(sizeof(double) *N);
	
	 
	
	rcode = 0;
	*niter = 0;
	siter = MAXITER;
	eps2 = sqrt(eps);
	gfdcode = 0;

	
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	fx = FUNCPT_EVAL(funcpt, xi, N);
	if (fx >= DBL_MAX || fx <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		rcode = 15;
	}
	if (fx != fx) {
		printf("Program Exiting as the function returns NaN");
		rcode = 15;
	}

	gfdcode = grad_fd(funcpt,funcgrad,xi,N,dx,eps2,jac);
	if (gfdcode == 15) {
		rcode = 15;
	}


	

	if (maxstep <= 0.0) {
		maxstep = 1000.0;
		dt1 = dt2 = 0.0;
		for (i = 0; i < N; ++i) {
			dt1 += dx[i] * dx[i];
			dt2 += dx[i] * xi[i] * dx[i] * xi[i];
		}

		dt1 = sqrt(dt1);
		dt2 = sqrt(dt2);

		if (dt1 > dt2) {
			maxstep *= dt1;
		}
		else {
			maxstep *= dt2;
		}
	}
	
	
	if (fabs(fx) > fabs(fsval)) {
			den = fabs(fx);
	} else {
			den = fabs(fsval);
	}
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xi[i]);
		} else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}
	
	stop0 = array_max_abs(scheck,N);
	
	if (stop0 <= gtol * 1e-03) {
		rcode = 1;
		for(i = 0; i < N;++i) {
			xf[i] = xi[i];
		}
	}
	
	
	inithess_lower(L,N,fx,fsval,dx);
	
	for(i = 0; i < N;++i) {
		xc[i] = xi[i];
	}
	
	while (rcode == 0 && *niter < siter) {
		*niter = *niter + 1;
		scale(jac,1,N,-1.0);
		
		linsolve_lower(L,N,jac,step);
		
		scale(jac,1,N,-1.0);

		retval = lnsrch(funcpt,xc,jac,step,N,dx,maxstep,stol,xf);

		fxf = FUNCPT_EVAL(funcpt, xf, N);
		if (fxf >= DBL_MAX || fxf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			rcode = 15;
			break;
		}
		if (fxf != fxf) {
			printf("Program Exiting as the function returns NaN");
			rcode = 15;
			break;
		}

		gfdcode = grad_fd(funcpt,funcgrad,xf,N,dx,eps2,jacf);
		if (gfdcode == 15) {
			rcode = 15;
			break;
		}
		rcode = stopcheck(fxf,N,xc,xf,jacf,dx,fsval,gtol,stol,retval);
		
		
		bfgs_factored(L,N,eps,xc,xf,jac,jacf);
		for(i = 0; i < N;++i) {
			xc[i] = xf[i];
			jac[i] = jacf[i];
		}
	}
	
	if (rcode == 0 && *niter >= siter) {
		rcode = 4;
	}
	
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	
	free(jac);
	free(hess);
	free(scheck);
	free(xc);
	free(L);
	free(step);
	free(jacf);
	return rcode;
}


/*
 * fminunc - 无约束优化统一入口 (SARIMA 使用 BFGS)
 */
int fminunc(custom_function *funcpt, custom_gradient *funcgrad, int N, double *xi,double maxstep, int method,double *xf) {
	int i,retval,MAXITER,niter,m;
	double fsval,eps,gtol,stol,ftol,xtol,delta;
	double *dx;

	dx = (double*) malloc(sizeof(double) * N);
		
	
	fsval = 1.0;
	MAXITER = 200*N;
	niter = 0;
	delta = -1.0; 
	eps = macheps(); 

	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}

	(void)method;
	(void)delta;
	(void)m;
	(void)ftol;
	(void)xtol;
	gtol = pow(eps, 1.0 / 3.0);
	stol = gtol * gtol;
	if (MAXITER < 1000)
		MAXITER = 1000;
	retval = bfgs_min(funcpt, funcgrad, xi, N, dx, fsval, maxstep, MAXITER, &niter, eps, gtol, stol, xf);



	free(dx);
	return retval;
}
