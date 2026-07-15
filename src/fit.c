// SPDX-License-Identifier: BSD-3-Clause
/*
 * fit.c - 季节 ARIMA 拟合 (emle)
 *
 * 功能: AS154/AS182 似然、CSS、STARMA/Kalman、MLE 目标函数与拟合例程。
 */

#include "fit.h"
#include "ts.h"
#include "optim.h"
#include "linalg.h"
#include "cstats.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>






/*
 * checkroots - 检查 AR/MA 平稳性与可逆性
 */
void checkroots(double *phi, int *p, double *theta, int *q, double *PHI, int *P, double *THETA, int *Q) {
	int ret;

	if (*p > 0) {
		ret = archeck(*p,phi);
		if (!ret) {
			printf("\nnon-stationary AR part\n");
			exit(-1);
		}
	}

	if (*q > 0) {
		invertroot(*q,theta);
	}

	if (*P > 0) {
		ret = archeck(*P, PHI);
		if (!ret) {
			printf("\nnon-stationary seasonal AR part\n");
			exit(-1);
		}
	}

	if (*Q > 0) {
		invertroot(*Q, THETA);
	}

}






/*
 * checkroots_cerr - checkroots 的 C 错误码封装
 */
static int checkroots_cerr(double *phi, int *p, double *theta, int *q, double *PHI, int *P, double *THETA, int *Q) {
	int ret,out;

	out = 1;

	if (*p > 0) {
		ret = archeck(*p,phi);
		if (!ret) {
			out = 10;
			printf("\nnon-stationary AR part\n");
			return out;
		}
	}

	if (*q > 0) {
		invertroot(*q,theta);
	}

	if (*P > 0) {
		ret = archeck(*P, PHI);
		if (!ret) {
			out = 12;
			printf("\nnon-stationary seasonal AR part\n");
			return out;
		}
	}

	if (*Q > 0) {
		invertroot(*Q, THETA);
	}

	return out;

}


/*
 * alik_css_seas_init - 初始化季节 CSS 似然工作区
 */
alik_css_seas_object alik_css_seas_init(int p, int d, int q, int s, int P, int D, int Q, int N)
{
	alik_css_seas_object obj = NULL;
	int i, r, t;

	r = p + s*P;

	t = q + s*Q + 1;

	if (r < t) {
		r = t;
	}
	else {
		t = r;
	}

	obj = (alik_css_seas_object)malloc(sizeof(struct alik_css_seas_set) + sizeof(double)* 2 * r + sizeof(double)* 3 * N);

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->s = s;
	obj->P = P;
	obj->D = D;
	obj->Q = Q;
	obj->N = N;
	obj->length = N;
	obj->pq = p + q + P + Q;

	obj->M = 0;

	if (d == 0 && D == 0) {
		obj->pq = obj->pq + 1;
		obj->M = 1;
	}

	obj->r = p + s * P;

	t = 1 + q + s*Q;
	if (obj->r < t) {
		obj->r = t;
	}
	t = obj->r;
	obj->offset = 2 * obj->r;

	for (i = 0; i < t; ++i) {
		obj->x[i] = 0.0;
	}

	for (i = 0; i < t; ++i) {
		obj->x[i + r] = 0.0;
	}
	obj->eps = macheps();
	obj->mean = 0.0;

	return obj;
}






/*
 * alik_seas_init - 初始化季节 MLE 似然工作区
 */
alik_seas_object alik_seas_init(int p, int d, int q, int s, int P, int D, int Q, int N) {
	alik_seas_object obj = NULL;
	int i,r, t;

	r = p + s*P;

	t = q + s*Q + 1;

	if (r < t) {
		r = t;
	}
	else {
		t = r;
	}

	obj = (alik_seas_object)malloc(sizeof(struct alik_seas_set) + sizeof(double) * 2 * r + sizeof(double)* 3 * N);

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->s = s;
	obj->P = P;
	obj->D = D;
	obj->Q = Q;
	obj->N = N;
	obj->length = N;
	obj->pq = p + q + P + Q;

	obj->M = 0;

	if (d == 0 && D == 0) {
		obj->pq = obj->pq + 1;
		obj->M = 1;
	}

	obj->r = p + s * P;

	t = 1 + q + s*Q;
	if (obj->r < t) {
		obj->r = t;
	}
	t = obj->r;
	obj->offset = 2 * obj->r;

	for (i = 0; i < t; ++i) {
		obj->x[i] = 0.0;
	}

	for (i = 0; i < t; ++i) {
		obj->x[i+r] = 0.0;
	}
	obj->eps = macheps();
	obj->mean = 0.0;

	return obj;
}


/*
 * inclu2 - Givens 旋转更新 Rbar (回归辅助)
 */
static int inclu2(int np, int nrbar, double weight, double *xnext, double *xrow, double ynext,
	double *d, double *rbar, double *thetab, double *ssqerr, double *recres, int *irank)
{
	int ifault, i, ithisr, i1, k;
	double y, wt, xi, di, dpi, cbar, sbar, xk, rbthis;

	(void)nrbar;
	y = ynext;
	wt = weight;

	for (i = 0; i < np; ++i) {
		xrow[i] = xnext[i];
	}

	*recres = 0.0;
	ifault = 1;
	if (wt <= 0.0) {
		return ifault;
	}
	ifault = 0;
	ithisr = 0;

	for (i = 0; i < np; ++i) {
		if (xrow[i] != 0.0) {
			xi = xrow[i];
			di = d[i];
			dpi = di + wt*xi*xi;
			d[i] = dpi;
			cbar = di / dpi;
			sbar = wt * xi / dpi;
			wt = cbar * wt;
			if (i != (np - 1)) {
				i1 = i + 1;
				for (k = i1; k < np; ++k) {
					xk = xrow[k];
					rbthis = rbar[ithisr];
					xrow[k] = xk - xi * rbthis;
					rbar[ithisr] = cbar * rbthis + sbar * xk;
					ithisr++;
				}
			}
			xk = y;
			y = xk - xi *thetab[i];
			thetab[i] = cbar * thetab[i] + sbar * xk;
			if (di == 0.0) {
				*irank = *irank + 1;
				return ifault;
			}
		}
		else {
			ithisr += (np - i - 1);
		}
	}
	*ssqerr = *ssqerr + wt*y*y;
	*recres = y * sqrt(wt);

	return ifault;
}






/*
 * regres - 从 Rbar 回代回归系数 (STARMA 用)
 */
static void regres(int np,int nrbar,double *rbar, double *thetab, double *beta) {
	int ithisr,i,im,i1,jm,j;
	double b1;
	
	ithisr = nrbar-1;
	im = np-1;
	for (i = 0; i < np; ++i) {
		b1 = thetab[im];
		if (im != (np - 1)) {
			i1 = i - 1;
			jm = np - 1;
			for (j = 0; j <= i1; ++j) {
				b1 = b1 - rbar[ithisr] * beta[jm];
				ithisr--;
				jm--;
			}
		}
		beta[im] = b1;
		im--;
	}
}









/*
 * starma - STARMA 初始条件与方差设置
 */
int starma(int ip, int iq,double *phi,double *theta,double *A,double *P,double *V) {
	int ifault,ir,np,nrbar,i,ind,j,irank;
	int npr, npr1, ind1, ind2, indj,indi,indn;
	double vj,ssqerr,phij,phii,ynext,weight,recres;
	double *thetab, *xnext, *xrow, *rbar;
	ifault = 0;
	

	

	if (ip < 0) {
		ifault = 1;
	}
	if (iq < 0) {
		ifault = 2;
	}

	if (ip*ip + iq*iq == 0) {
		ifault = 4;
	}

	ir = iq + 1;

	if (ir < ip) {
		ir = ip;
	}

	np = (ir * (ir + 1)) / 2;
	nrbar = (np * (np - 1)) / 2;

	if (ir == 1) {
		ifault = 8;
	}
	
	if (ifault != 0) {
		return ifault;
	}

	xnext = (double*)malloc(sizeof(double)* np);
	thetab = (double*)malloc(sizeof(double)* np);
	xrow = (double*)malloc(sizeof(double)* np);
	rbar = (double*)malloc(sizeof(double)* nrbar);

	for (i = 1; i < ir;++i) {
		A[i] = 0.0;
		if (i >= ip) {
			phi[i] = 0.0;
		}
		V[i] = 0.0;
		if (i < iq+1) {
			V[i] = theta[i - 1];
		}
	}

	A[0] = 0.0;
	if (ip == 0) {
		phi[0] = 0.0;
	}
	V[0] = 1.0;
	ind = ir;

	for (j = 1; j < ir; ++j) {
		vj = V[j];
		for (i = j; i < ir; ++i) {
			V[ind] = V[i] * vj;
			ind++;
		}
	}

	if (ip != 0) {
		
		irank = 0;
		ssqerr = 0.0;

		for (i = 0; i < nrbar; ++i) {
			rbar[i] = 0.0;
		}
		for (i = 0; i < np; ++i) {
			P[i] = 0.0;
			thetab[i] = 0.0;
			xnext[i] = 0.0;
		}
		ind = 0;
		ind1 = -1;
		npr = np - ir;
		npr1 = npr + 1;
		indj = npr;
		ind2 = npr-1;

		for (j = 0; j < ir; ++j) {
			phij = phi[j];
			xnext[indj] = 0.0;
			indj = indj + 1;
			indi = npr1 + j;

			for (i = j; i < ir; ++i) {
				ynext = V[ind];
				ind++;
				phii = phi[i];
				if (j != (ir - 1)) {
					xnext[indj] = -phii;
					if (i != (ir - 1)) {
						xnext[indi] -= phij;
						ind1 ++;
						xnext[ind1] = -1.0;
					}
				}
				xnext[npr] = -phii*phij;
				ind2++;
				if (ind2 >= np) {
					ind2 = 0;
				}
				xnext[ind2] += 1.0;
				weight = 1.0;
				(void)inclu2(np, nrbar, weight, xnext, xrow, ynext, P, rbar, thetab, &ssqerr, &recres, &irank);
				
				xnext[ind2] = 0.0;
				if (i != (ir - 1)) {
					xnext[indi] = 0.0;
					indi++;
					xnext[ind1] = 0.0;
				}
			}
		}
		
		regres(np, nrbar, rbar, thetab, P);
		
		ind = npr;
		for (i = 0; i < ir; ++i) {
			ind++;
			xnext[i] = P[ind - 1];
		}
		ind = np;
		ind1 = npr;

		for (i = 0; i < npr; ++i) {
			P[ind - 1] = P[ind1 - 1];
			ind--;
			ind1--;
		}

		for (i = 0; i < ir; ++i) {
			P[i] = xnext[i];
		}
		
	}
	else {

		

		indn = np;
		ind = np;
		for (i = 0; i < ir; i++) {
			for (j = 0; j <= i; j++) {
				ind--;
				P[ind] = V[ind];
				if (j != 0)  {
					indn--;
					P[ind] += P[indn];
				}
			}
		}
	}

	free(xnext);
	free(thetab);
	free(xrow);
	free(rbar);
	return ifault;
}







/*
 * karma - Kalman 滤波递推一步
 */
void karma(int ip,int iq,double *phi,double *theta,double *A,double *P,double*V,int N,
	double *W,double *resid,double *sumlog,double *ssq,int iupd,double delta,int *iter,int *nit) {
	int ir,i,j,ir1,inde,swtch,ind,indn,l,ii,indw;
	double wnext,A1,ft,ut,g,et;
	double *E;

	(void)delta;
	ir = iq + 1;
	swtch = 0;
	*iter = 0;

	if (ir < ip) {
		ir = ip;
	}
	E = (double*)malloc(sizeof(double)* ir);

	ir1 = ir - 1;
	for (i = 0; i < ir; ++i) {
		E[i] = 0.0;
	}
	inde = 0;

	if (*nit == 0) {
		for (i = 0; i < N; ++i) {
			wnext = W[i];
			if (iupd != 1 || i > 0) {
				
				A1 = A[0];
				if (ir != 1) {
					for (j = 0; j < ir1; ++j) {
						A[j] = A[j + 1];
					}
				}
				A[ir - 1] = 0.0;
				if (ip != 0) {
					for (j = 0; j < ip; ++j) {
						A[j] += phi[j] * A1;
					}
				}
				ind = -1;
				indn = ir-1;
				for (l = 0; l < ir; ++l) {
					for (j = l; j < ir; ++j) {
						ind++;
						P[ind] = V[ind];
						if (j != (ir - 1)) {
							indn++;
							P[ind] += P[indn];
						}
					}
				}
			}
			ft = fabs(P[0]);
			ut = wnext - A[0];
			if (ir != 1) {
				ind = ir;
				for (j = 1; j < ir; ++j) {
					g = P[j] / ft;
					A[j] += g*ut;
					for (l = j; l < ir; ++l) {
						P[ind] -= g*P[l];
						ind++;
					}
				}
			}
			A[0] = wnext;
			for (l = 0; l < ir; ++l) {
				P[l] = 0.0;
			}
			resid[i] = ut / sqrt(ft);
			E[inde] = resid[i];
			inde++;
			if (inde >= iq) {
				inde = 0;
			}
			*ssq += (ut*ut) / ft;
			*sumlog += log(ft);
			*iter = *iter + 1;;
			


		}
		if (swtch == 0) {
			*nit = N;
		}
	}
	else {
		i = 0;
		*nit = i;
		for (ii = i; ii < N; ++ii) {
			et = W[ii];
			indw = ii;
			if (ip != 0) {
				for (j = 0; j < ip; ++j) {
					indw--;
					if (indw >= 0) {
						et -= phi[j] * W[indw];
					}
				}
			}
			if (iq != 0) {
				for (j = 0; j < iq; ++i) {
					inde--;
					if (inde == -1) {
						inde = iq-1;
					}
					et -= theta[j] * E[inde];
				}
			}
			E[inde] = et;
			resid[ii] = et;
			*ssq += et*et;
			*iter = *iter + 1;
			inde++;
			if (inde >= iq) {
				inde = 0;
			}
		}

	}
	if (*iter == 0) {
		*iter = 1;
	}
	free(E);
}









/*
 * forkal - Kalman 型多步预测与 AMSE
 */
int forkal(int ip,int iq,int id,double *phi,double*theta,double *delta,int N,double *W,double *resid,int il,double *Y,double *AMSE) {
	int ifault,ir,np,k,ird,irz;
	double *A,*P,*V,*store,*xrow;
	double zero, one, two,AA,del,sumlog,ssq,sigma,A1,dt,phij,phijdt,phii,AMS;
	int i,  j, ll, nt, nj, idk, iid,iupd,nit,iter,ind,irj;
	int ir2, ir1, id2r, id2r1, id1,idd1,idd2,i45,idrr1,iddr,jkl,jkl1,id2r2,ibc,l,iri1,jj;
	int lk, lk1,jklj,iri,kk1,k1,kk,kkk,ind1,ind2,jrj,j1,jrk;
	
	ir = iq + 1;

	if (ir < ip) {
		ir = ip;
	}
	np = (ir * (ir + 1)) / 2;
	ird = ir + id;
	irz = (ird * (ird + 1)) / 2;
	zero = 0.0;
	one = 1.0;
	two = 2.0;

	ifault = 0;

	if (ip < 0) {
		ifault = 1;
	}
	if (iq < 0) {
		ifault = ifault + 2;
	}
	if (ip*ip + iq*iq == 0) {
		ifault = 4;
	}
	if (id < 0) {
		ifault = 8;
	}
	if (il < 1) {
		ifault = 11;
	}
	if (ifault != 0) {
		return ifault;
	}

	A = (double*)malloc(sizeof(double)* ird);
	P = (double*)malloc(sizeof(double)* irz);
	V = (double*)malloc(sizeof(double)* np); 
	store = (double*)malloc(sizeof(double)* ird);
	xrow = (double*)malloc(sizeof(double)* np);

	

	for (i = 0; i < ird; ++i) {
		A[i] = zero;
		store[i] = zero;
	}
	for (i = 0; i < irz; ++i) {
		P[i] = zero;
	}
	for (i = 0; i < np; ++i) {
		V[i] = zero;
		xrow[i] = zero;
	}
	A[0] = zero;
	V[0] = one;
	
	

	if (ir == 1) {
		*P = 1.0 / (1.0 - phi[0] * phi[0]);
	}
	else {
		starma(ip, iq, phi, theta,A, P, V);
	}

	

	nt = N - id;

	if (id != 0) {
		for (j = 1; j <= id; ++j) {
			nj = N - j;
			store[j-1] = W[nj-1];
		}
		for (i = 1; i <= nt; ++i) {
			AA = zero;
			for (k = 1; k <= id; ++k) {
				idk = id + i - k;
				AA -= delta[k - 1] * W[idk - 1];
			}
			iid = i + id;
			W[i - 1] = W[iid - 1] + AA;
		}
	}

	
	sumlog = ssq = zero;
	del = -1.0;
	iupd = 1;
	nit = 0;
	iter = 0;

	karma(ip, iq, phi, theta, A, P, V, nt, W, resid, &sumlog, &ssq, iupd, del,&iter, &nit);

	

	sigma = zero;

	for (j = 0; j < nt; ++j) {
		sigma += resid[j] * resid[j];
	}
	sigma = sigma / nt;
	

	

	if (id != 0) {
		for (i = 1; i <= np; ++i) {
			xrow[i-1] = P[i-1];
		}
		for (i = 1; i <= irz; ++i) {
			P[i-1] = zero;
		}
		ind = 0;
		for (j = 1; j <= ir; ++j) {
			k = (j - 1) * (id + ir + 1) - (j - 1) * j / 2;
			for (i = j; i <= ir; ++i) {
				ind++;
				k++;
				P[k - 1] = xrow[ind-1];
			}
		}

		for (j = 1; j <= id; ++j) {
			irj = ir + j;
			A[irj-1] = store[j-1];
		}

	}

	
	ir2 = ir + 1;
	ir1 = ir - 1;
	id1 = id - 1;
	id2r = 2 * ird;
	id2r1 = id2r - 1;
	idd1 = 2 * id + 1;
	idd2 = idd1 + 1;
	i45 = id2r + 1;
	idrr1 = ird + 1;
	iddr = 2 * id + ir;
	jkl = ir * (iddr + 1) / 2;
	jkl1 = jkl + 1;
	id2r2 = id2r + 2;
	ibc = ir * (i45 - ir) / 2;

	for (l = 1; l <= il; ++l) {
		
		A1 = A[0];
		if (ir != 1) {
			for (i = 1; i <= ir1; ++i) {
				A[i-1] = A[i];
			}
		}
		A[ir-1] = zero;
		if (ip != 0) {
			for (j = 1; j <= ip; ++j) {
				A[j-1] += phi[j-1] * A1;
			}
		}
		if (id != 0) {
			for (j = 1; j <= id; ++j) {
				irj = ir + j;
				A1 += delta[j-1] * A[irj-1];
			}
			if (id >= 2) {
				for (i = 1; i <= id1; ++i) {
					iri1 = ird - i;
					A[iri1] = A[iri1-1];
				}
			}
			A[ir2 - 1] = A1;
		}

		

		

		if (id != 0) {
			for (i = 1; i <= id; ++i) {
				store[i-1] = zero;
				for (j = 1; j <= id; ++j) {
					ll = imax(i, j);
					k = imin(i, j);
					jj = jkl + (ll - k) + 1 + (k - 1) * (idd2 - k) / 2;
					store[i-1] += delta[j-1] * P[jj-1];
				}
			}

			if (id != 1) {
				for (j = 1; j <= id1; ++j) {
					jj = id - j;
					lk = (jj - 1) * (idd2 - jj) / 2 + jkl;
					lk1 = jj * (idd1 - jj) / 2 + jkl;
					for (i = 1; i <= j; ++i) {
						lk = lk + 1;
						lk1 = lk1 + 1;
						P[lk1-1] = P[lk-1];
					}
				}
				for (j = 1; j <= id1; ++j) {
					jklj = jkl1 + j;
					irj = ir + j;
					P[jklj - 1] = store[j - 1] + P[irj-1];
				}
			}

			P[jkl1 - 1] = P[0];

			for (i = 1; i <= id; ++i) {
				iri = ir + i;
				P[jkl1 - 1] += delta[i - 1] * (store[i - 1] + two * P[iri - 1]);
			}

			for (i = 1; i <= id; ++i) {
				iri = ir + i;
				store[i - 1] = P[iri - 1];
			}

			for (j = 1; j <= ir; ++j) {
				kk1 = j * (id2r1 - j) / 2 + ir;
				k1 = (j - 1) * (id2r - j) / 2 + ir;
				for (i = 1; i <= id; ++i) {
					kk = kk1 + i;
					k = k1 + i;
					P[k - 1] = phi[j - 1] * store[i - 1];
					if (j != ir) {
						P[k - 1] += P[kk - 1];
					}
				}
			}

			for (j = 1; j <= ir; ++j) {
				store[j - 1] = zero;
				kkk = j * (i45 - j) / 2 - id;
				for (i = 1; i <= id; ++i) {
					kkk++;
					store[j - 1] += delta[i - 1] * P[kkk - 1];
				}
			}

			if (id != 1) {
				for (j = 1; j <= ir; ++j) {
					k = j * idrr1 - j * (j + 1) / 2 + 1;
					for (i = 1; i <= id1; ++i) {
						k--;
						P[k - 1] = P[k - 2];
					}
				}
			}

			for (j = 1; j <= ir; ++j) {
				k = (j - 1) * (id2r - j) / 2 + ir + 1;
				P[k - 1] = store[j - 1] + phi[j - 1] * P[0];
				if (j < ir) {
					P[k - 1] += P[j];
				}
			}

		}

		for (i = 0; i < ir; ++i) {
			store[i] = P[i];
		}

		ind = 0;
		dt = P[0];
		for (j = 1; j <= ir; ++j) {
			phij = phi[j - 1];
			phijdt = phij * dt;
			ind2 = (j - 1) * (id2r2 - j) / 2;
			ind1 = j * (i45 - j) / 2;
			for (i = j; i <= ir; ++i) {
				ind++;
				ind2++;
				phii = phi[i - 1];
				P[ind2 - 1] = V[ind - 1] + phii * phijdt;
				if (j < ir) {
					P[ind2 - 1] += store[j] * phii;
				}
				if (i != ir) {
					ind1++;
					P[ind2 - 1] += store[i] * phij + P[ind1 - 1];
				}
			}
		}

		

		Y[l - 1] = A[0];
		if (id != 0) {
			for (j = 1; j <= id; ++j) {
				irj = ir + j;
				Y[l - 1] += A[irj - 1] * delta[j - 1];
			}
			
		}
		AMS = P[0];
		if (id != 0) {
			for (j = 1; j <= id; ++j) {
				jrj = ibc + (j - 1) * (idd2 - j) / 2;
				irj = ir + j;
				AMS += (two * delta[j - 1] * P[irj - 1] + P[jrj] * delta[j - 1] * delta[j - 1]);
			}
			if (id != 1) {
				for (j = 1; j <= id1; ++j) {
					j1 = j + 1;
					jrk = ibc + 1 + (j - 1) * (idd2 - j) / 2;
					for (i = j1; i <= id; ++i) {
						jrk++;
						AMS += two * delta[i-1] * delta[j-1] * P[jrk-1];
					}
				}
			}
		}
		AMSE[l - 1] = AMS * sigma;

	}
	
	free(A);
	free(P);
	free(V);
	free(store);
	free(xrow);
	return ifault;
}






/*
 * fcss_seas - CSS 目标函数 (季节)
 */
double fcss_seas(double *b, int pq, void *params) {
	double value, ssq, temp;
	int p, q, ps, qs, s, offset, jm;
	int ip, iq, i, j, N, iter, ncond;
	double *phi, *theta;

	(void)pq;
	value = ssq = 0.0;
	
	alik_css_seas_object obj = (alik_css_seas_object)params;

	ip = obj->p + (obj->s * obj->P);
	iq = obj->q + (obj->s * obj->Q);

	N = obj->N;
	p = obj->p;
	ps = obj->P;
	q = obj->q;
	qs = obj->Q;
	s = obj->s;
	offset = obj->offset;
	ncond = p + s* ps;

	

	phi = (double*)malloc(sizeof(double)* ip);
	theta = (double*)malloc(sizeof(double)* iq);



	for (i = 0; i < p; ++i) {
		phi[i] = b[i];
	}

	for (i = p; i < ip; ++i) {
		phi[i] = 0.0;
	}

	for (i = 0; i < q; ++i) {
		theta[i] = b[i + p];
	}

	for (i = q; i < iq; ++i) {
		theta[i] = 0.0;
	}

	for (j = 0; j < ps; ++j) {
		phi[(j + 1)*s - 1] += b[p + q + j];
		for (i = 0; i < p; ++i) {
			phi[(j + 1)*s + i] -= b[i] * b[p + q + j];
		}
	}

	for (j = 0; j < qs; ++j) {
		theta[(j + 1)*s - 1] += b[p + q + ps + j];
		for (i = 0; i < q; ++i) {
			theta[(j + 1)*s + i] += b[i + p] * b[p + q + ps + j];
		}
	}
	
	
	for (i = 0; i < ncond; ++i) {
		obj->x[offset+N + i] = 0.0;
	}

	if (obj->M == 1) {
		for (i = 0; i < N; ++i) {
			obj->x[offset+i] = obj->x[offset+ 2 * N + i] - b[p+q+ps+qs];
		}
	}
	iter = 0;
	for (i = ncond; i < N; ++i) {
		iter++;
		temp = obj->x[offset+i];

		for (j = 0; j < ip; ++j) {
			temp = temp - phi[j] * obj->x[offset + i - j - 1];
		}

		if (i - ncond < iq) {
			jm = i - ncond;
		}
		else {
			jm = iq;
		}

		for (j = 0; j < jm; ++j) {
			temp = temp - theta[j] * obj->x[offset+N + i - j - 1];
		}

		obj->x[offset+N + i] = temp;
		ssq += temp * temp;
	}
	obj->ssq = ssq;
	value = 0.5 * log(ssq / (double)iter);
	obj->loglik = value;

	

	free(phi);
	free(theta);
	return value;
}






/*
 * css_seas - 季节 ARIMA 条件平方和 CSS
 */
int css_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,
	double *phi, double *theta, double *PHI, double *THETA, double *wmean, double *var,double *loglik,double *hess) {
	int i, pq, retval, offset,ret;
	double *b, *tf, *x, *inp2,*dx,*thess;
	int *ipiv;
	double maxstep;
	alik_css_seas_object obj;
	

	inp2 = (double*)malloc(sizeof(double)* (N - s*D));

	maxstep = 1.0;

	obj = alik_css_seas_init(p, d, q, s, P, D, Q, N);

	pq = obj->pq;
	b = (double*)malloc(sizeof(double)* pq);
	tf = (double*)malloc(sizeof(double)* pq);
	thess = (double*)malloc(sizeof(double)* pq*pq);
	dx = (double*)malloc(sizeof(double)* pq);
	ipiv = (int*)malloc(sizeof(int)* pq);
	

	if (D > 0) {
		N = diffs(inp, N, D, s, inp2);
	}
	else {
		for (i = 0; i < N; ++i) {
			inp2[i] = inp[i];
		}
	}

	x = (double*)malloc(sizeof(double)* (N - d));

	if (d > 0) {
		N = diff(inp2, N, d, x); 
	}
	else {
		for (i = 0; i < N; ++i) {
			x[i] = inp2[i];
		}
	}

	obj->N = N;

	offset = obj->offset;
	

	for (i = 0; i < p; ++i) {
		b[i] = 0.0;
	}
	for (i = 0; i < q; ++i) {
		b[p + i] = 0.0;
	}
	for (i = 0; i < P; ++i) {
		b[p + q + i] = 0.0;
	}
	for (i = 0; i < Q; ++i) {
		b[p + q + P + i] = 0.0;
	}

	if (obj->M == 1) {
		b[p + q + P + Q] = mean(x, N);
		*wmean = b[p + q + P + Q];
	}
	else {
		*wmean = 0.0;
	}

	obj->mean = *wmean;

	

	for (i = 0; i < N; ++i) {
		obj->x[offset + i] = obj->x[offset + 2 * N + i] = x[i];
	}
	for (i = N; i < 2 * N; ++i) {
		obj->x[offset + i] = 0.0;
	}

	custom_function css_min = { fcss_seas, obj };
	retval = fminunc(&css_min, NULL, pq, b, maxstep, optmethod, tf);
	if (retval == 0) {
		ret = 0;
	}
	else if (retval == 15) {
		ret = 15;
	}
	else if (retval == 4) {
		ret = 4;
	}
	else {
		ret = 1;
	}

	for (i = 0; i < pq; ++i) {
		dx[i] = 1.0;
	}

	hessian_fd(&css_min, tf, pq, dx, obj->eps, hess);

	mtranspose(hess, pq, pq, thess);

	for (i = 0; i < pq*pq; ++i) {
		thess[i] = (N - d - s*D) * 0.5 * (hess[i] + thess[i]);
	}


	ludecomp(thess, pq, ipiv);
	minverse(thess, pq, ipiv, hess);

	for (i = 0; i < p; ++i) {
		phi[i] = tf[i];
	}
	for (i = 0; i < q; ++i) {
		theta[i] = -tf[p + i];
	}
	for (i = 0; i < P; ++i) {
		PHI[i] = tf[p + q + i];
	}
	for (i = 0; i < Q; ++i) {
		THETA[i] = -tf[p + q + P + i];
	}

	if (obj->M == 1) {
		*wmean = tf[p + q + Q + P];
	}
	else {
		*wmean = 0.0;
	}
	
	

	*var = (obj->ssq) / (double)N;
	*loglik = obj->loglik;
	
	

	free(b);
	free(tf);
	free(inp2);
	free(x);
	free(thess);
	free(dx);
	free(ipiv);
	free_alik_css_seas(obj);
	return ret;
}






/*
 * fas154_seas - AS154 MLE 目标函数 (季节)
 */
double fas154_seas(double *b, int pq, void *params) {
	double value, ssq, sumlog, delta;
	int p, q, ps, qs,s,offset;
	int ip, iq, ir, i,j, np, N, iupd, nit, iter;
	double *phi, *theta, *A, *P, *V;

	(void)pq;

	value = ssq = sumlog = 0.0;
	alik_seas_object obj = (alik_seas_object)params;

	ip = obj->p + (obj->s * obj->P);
	iq = obj->q + (obj->s * obj->Q);
	ir = obj->r;
	N = obj->N;
	p = obj->p;
	ps = obj->P;
	q = obj->q;
	qs = obj->Q;
	s = obj->s;
	offset = obj->offset;
	np = (ir * (ir + 1)) / 2;

	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	A = (double*)malloc(sizeof(double)* ir);
	P = (double*)malloc(sizeof(double)* np);
	V = (double*)malloc(sizeof(double)* np);

	for (i = 0; i < ir; ++i) {
		phi[i] = 0.0;
		theta[i] = 0.0;
	}

	for (i = 0; i < p; ++i) {
		phi[i] = b[i];
	}

	for (i = 0; i < q; ++i) {
		theta[i] = b[i + p];
	}
	

	for (j = 0; j < ps; ++j) {
		phi[(j + 1)*s - 1] += b[p + q + j];
		for (i = 0; i < p; ++i) {
			phi[(j + 1)*s + i] -= b[i] * b[p + q + j];
		}
	}

	for (j = 0; j < qs; ++j) {
		theta[(j + 1)*s - 1] += b[p + q + ps + j];
		for (i = 0; i < q; ++i) {
			theta[(j + 1)*s + i] += b[i + p] * b[p + q + ps + j];
		}
	}
	if (obj->M == 1) {
		for (i = 0; i < N; ++i) {
			obj->x[offset + i] = obj->x[offset + 2 * N + i] - b[p + q + ps + qs];
		}
	}
	
	
	
	iupd = 0;
	

	if (ip == 1 && iq == 0) {
		*V = 1.0;
		*A = 0.0;
		*P = 1.0 / (1.0 - phi[0] * phi[0]);
	}
	else {
		iupd = 1;
		(void)starma(ip, iq, phi, theta, A, P, V);
	}

	nit = 0;
	delta = 0.001;
	karma(ip, iq, phi, theta, A, P, V, N, obj->x + offset, obj->x + offset + N, &sumlog, &ssq, iupd, delta, &iter, &nit);
	obj->ssq = ssq;
	

	value = 0.5 * (sumlog / (double)iter + log(ssq / (double)iter));
	obj->loglik = value;
	

	free(phi);
	free(theta);
	free(A);
	free(P);
	free(V);
	return value;
}






/*
 * as154_seas - 季节 ARIMA AS154 精确 MLE
 */
int as154_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,double *phi, double *theta, 
	double *PHI, double *THETA, double *wmean,double *var,double *loglik,double *hess,int cssml) {
	int i, pq, retval, length, offset, ret, nd, ERR;
	double *b, *tf, *x, *inp2, *dx, *thess;
	int *ipiv;
	double maxstep, coeff, sigma;
	alik_seas_object obj;

	inp2 = (double*)malloc(sizeof(double)* (N - s*D));

	length = N;
	maxstep = 1.0;
	coeff = 0.0;
	sigma = 1.0;
	(void)sigma;

	if (cssml == 1) {
		css_seas(inp, N, optmethod, p, d, q, s, P, D, Q, phi, theta, PHI, THETA, wmean, var, loglik, hess);

		ERR = checkroots_cerr(phi, &p, theta, &q, PHI, &P, THETA, &Q);
		if (ERR == 10 || ERR == 12) return ERR;

		coeff = *wmean;
	}
	else {
		for (i = 0; i < p; ++i) phi[i] = 0.0;
		for (i = 0; i < q; ++i) theta[i] = 0.0;
		for (i = 0; i < P; ++i) PHI[i] = 0.0;
		for (i = 0; i < Q; ++i) THETA[i] = 0.0;
	}


	

	nd = d + D;

	if (D > 0) {
		N = diffs(inp, N, D, s, inp2);
	}
	else {
		for (i = 0; i < N; ++i) {
			inp2[i] = inp[i];
		}
	}

	if (p + q + P + Q == 0 && d == 0 && D == 0 && cssml == 1) {
		free(inp2);
		return 1;
	}

	x = (double*)malloc(sizeof(double)* (N - d));

	if (d > 0) {
		N = diff(inp2, N, d, x); 
	}
	else {
		for (i = 0; i < N; ++i) {
			x[i] = inp2[i];
		}
	}

	if (nd == 0) {
		double sum = 0.0, ss = 0.0;
		for (i = 0; i < N; ++i)
			sum += x[i];
		coeff = sum / (double)N;
		for (i = 0; i < N; ++i) {
			double e = x[i] - coeff;
			ss += e * e;
		}
		sigma = (N > 1) ? 10.0 * sqrt(ss / (double)(N - 1) / (double)N) : 1.0;
	}

	obj = alik_seas_init(p, d, q, s, P, D, Q, N);
	pq = obj->pq;
	b = (double*)malloc(sizeof(double)* pq);
	tf = (double*)malloc(sizeof(double)* pq);
	thess = (double*)malloc(sizeof(double)* pq*pq);
	dx = (double*)malloc(sizeof(double)* pq);
	ipiv = (int*)malloc(sizeof(int)* pq);

	obj->N = N;

	offset = obj->offset;
	for (i = 0; i < p; ++i) {
		b[i] = phi[i];
	}
	for (i = 0; i < q; ++i) {
		b[p + i] = -theta[i];
	}
	for (i = 0; i < P; ++i) {
		b[p + q + i] = PHI[i];
	}
	for (i = 0; i < Q; ++i) {
		b[p + q + P + i] = -THETA[i];
	}

	if (obj->M == 1) {
		b[p + q + P + Q] = coeff;
	}

	obj->mean = coeff;

	

	for (i = 0; i < N; ++i) {
		obj->x[offset + i] = obj->x[offset + 2 * N + i] = x[i];
	}
	for (i = N; i < 2 * N; ++i) {
		obj->x[offset + i] = 0.0;
	}
	

	custom_function as154_min = { fas154_seas, obj };
	retval = fminunc(&as154_min, NULL, pq, b, maxstep, optmethod, tf);
	if (retval == 0) {
		ret = 0;
	}
	else if (retval == 15) {
		ret = 15;
	}
	else if (retval == 4) {
		ret = 4;
	}
	else {
		ret = 1;
	}

	for (i = 0; i < pq; ++i) {
		dx[i] = 1.0;
	}

	

	if (q > 0) {
		invertroot(q,tf+p);
	}

	if (Q > 0) {
		invertroot(Q,tf+p+q+P);
	}

	

	hessian_fd(&as154_min, tf, pq, dx, obj->eps, hess);

	mtranspose(hess, pq, pq, thess);

	for (i = 0; i < pq*pq; ++i) {
		thess[i] = (length - d - s*D) * 0.5 * (hess[i] + thess[i]);
	}

	

	ludecomp(thess, pq, ipiv);
	minverse(thess, pq, ipiv, hess);

	


	for (i = 0; i < p; ++i) {
		phi[i] = tf[i];
	}
	for (i = 0; i < q; ++i) {
		theta[i] = -tf[p + i];
	}
	for (i = 0; i < P; ++i) {
		PHI[i] = tf[p + q + i];
	}
	for (i = 0; i < Q; ++i) {
		THETA[i] = -tf[p + q + P + i];
	}

	if (obj->M == 1) {
		*wmean = tf[p + q + Q + P];
	}
	else {
		*wmean = 0.0;
	}

	*var = (obj->ssq) / (double) N;
	*loglik = obj->loglik;
	
	

	free(b);
	free(tf);
	free(inp2);
	free(x);
	free(dx);
	free(thess);
	free(ipiv);
	free_alik_seas(obj);
	return ret;
}






/*
 * free_alik_seas - 释放季节 MLE 工作区
 */
void free_alik_seas(alik_seas_object object)
{
	free(object);
}


/*
 * free_alik_css_seas - 释放季节 CSS 工作区
 */
void free_alik_css_seas(alik_css_seas_object object)
{
	free(object);
}
