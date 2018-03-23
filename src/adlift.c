#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <R_ext/RS.h> ///* for F77_CALL */

#include "adlift.h"

/*
	Bottlenecks:
		getridd
		myt
		linearpred_nn_not_one
		getridi
		getridd
*/

#define MULTIPLY2(a) a << 1
#define MULTIPLY3(a) a + (a << 1)
#define MULTIPLY4(a) a << 2

/* ********** */

void adaptneigh(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	 int *inter, int *nn, double *weights, int *scheme, int *clo, int *index,
	 int *neighbours, int *N, int *docoeff)
{
	int i, tot, min1, min2, ip;
	int minindex;
	double min, *absdetails, cr;

	double *w;
	int *nbrs1;
	int *nbrs2, *index2;
	int one = 1;

	min1 = ((*N - 1) <= *neighbours) ? (*N - 1) : *neighbours;
	min2 = ((*N - 1) <= MULTIPLY2(*neighbours)) ? (*N - 1) : MULTIPLY2(*neighbours);

	tot = min1 + min2;
	absdetails = calloc(tot, sizeof(double));
	cr = *(coeff + *remove - 1);

	*clo = 0;

	/*nbrs2=calloc(2**neighbours,sizeof(int));
	 index2=calloc(2**neighbours,sizeof(int));*/
	/* let's use temp variables after all, and then output the proper arrays in the pointers in the last call. */

	for (i = 0; i < min1; i++) {
		ip = i + 1;
		/*    nbrs2=calloc(2**neighbours,sizeof(int));
		 index2=calloc(2**neighbours,sizeof(int));*/
		nbrs2 = calloc(MULTIPLY2(ip), sizeof(int));
		index2 = calloc(MULTIPLY2(ip), sizeof(int));
		getnbrs(X, remove, pointsin, N, &ip, clo, nbrs2, index2, nn);

		nbrs1 = calloc(*nn, sizeof(int));

		mycpyi(nbrs2, nn, nbrs1);

		free(nbrs2);
		/* MAN: 2/4/09 */
		free(index2);
		w = calloc(*nn, sizeof(double));

		adaptpred(pointsin, X, coeff, nbrs1, remove, inter, nn, w, scheme,
				&one);

		*(absdetails + i) = fabs(*(coeff + *remove - 1));
		*(coeff + *remove - 1) = cr; /*  reset of coeff */
		/* others (nbrs,inter,nn,weights,sch) shouldn't matter or get reset anyway. */

		free(w);
		free(nbrs1);
	}

	*clo = 1;

	for (i = 0; i < min2; i++) {
		/*    nbrs2=calloc(*neighbours,sizeof(int));
		 index2=calloc(*neighbours,sizeof(int));*/
		ip = i + 1;
		nbrs2 = calloc(ip, sizeof(int));
		index2 = calloc(ip, sizeof(int));
		getnbrs(X, remove, pointsin, N, &ip, clo, nbrs2, index2, nn);

		nbrs1 = calloc(*nn, sizeof(int));

		mycpyi(nbrs2, nn, nbrs1);

		free(nbrs2);
		/* MAN: 2/4/09 */
		free(index2);

		w = calloc(*nn, sizeof(double));

		adaptpred(pointsin, X, coeff, nbrs1, remove, inter, nn, w, scheme,
				&one);

		*(absdetails + min1 + i) = fabs(*(coeff + *remove - 1));
		*(coeff + *remove - 1) = cr; /*  reset of coeff */
		/* others (nbrs,inter,nn,weights,sch) shouldn't matter or get reset anyway. */

		free(w);
		free(nbrs1);
	}

	mymind(absdetails, &tot, &min, &minindex);

	free(absdetails);

	*(coeff + *remove - 1) = cr; /* should be reset already,but just in case something wierd happens  */

	if (minindex <= min1) {
		*clo = 0;
		ip = minindex;
	} else {
		*clo = 1; /* should be already */
		ip = minindex - min1;
	}

	getnbrs(X, remove, pointsin, N, &ip, clo, nbrs, index, nn); /* should
	 output "true" nbrs,nn and index */
	adaptpred(pointsin, X, coeff, nbrs, remove, inter, nn, weights, scheme,
			docoeff); /* should get correct coeff,inter,weights and scheme */

}

/* ********* */

void adaptpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	 int *inter, int *nn, double *weights, int *scheme, int *docoeff)
{

	int one = 1, six = 6, minindex;
	double *tmpweights, cr, min;
	double *details;

	/*void linearpred();
	void quadpred();
	void cubicpred();
	//void mycpyd();
	void mymind();
	void myrbind();*/

	cr = *(coeff + *remove - 1);

	*inter = 0;

	details = calloc(6, sizeof(double));

	tmpweights = calloc(*nn, sizeof(double));

	linearpred(pointsin, X, coeff, nbrs, remove, inter, nn, tmpweights, &one);
	*details = fabs(*(coeff + *remove - 1));
	*(coeff + *remove - 1) = cr;

	quadpred(pointsin, X, coeff, nbrs, remove, inter, nn, tmpweights, &one);
	*(details + 1) = fabs(*(coeff + *remove - 1));
	*(coeff + *remove - 1) = cr;

	cubicpred(pointsin, X, coeff, nbrs, remove, inter, nn, tmpweights, &one);
	*(details + 2) = fabs(*(coeff + *remove - 1));
	*(coeff + *remove - 1) = cr;

	*inter = 1;

	linearpred(pointsin, X, coeff, nbrs, remove, inter, nn, tmpweights, &one);
	*(details + 3) = fabs(*(coeff + *remove - 1));
	*(coeff + *remove - 1) = cr;

	quadpred(pointsin, X, coeff, nbrs, remove, inter, nn, tmpweights, &one);
	*(details + 4) = fabs(*(coeff + *remove - 1));
	*(coeff + *remove - 1) = cr;

	cubicpred(pointsin, X, coeff, nbrs, remove, inter, nn, tmpweights, &one);
	*(details + 5) = fabs(*(coeff + *remove - 1));
	*(coeff + *remove - 1) = cr;

	mymind(details, &six, &min, &minindex);

	free(details);
	free(tmpweights);

	*(coeff + *remove - 1) = cr;

	*inter = 1;
	if (minindex <= 3) {
		*inter = 0;
	}

	switch (minindex) {
	case 1:
		*scheme = 1;
		break;

	case 2:
		*scheme = 2;
		break;

	case 3:
		*scheme = 3;
		break;

	case 4:
		*scheme = 1;
		break;

	case 5:
		*scheme = 2;
		break;

	case 6:
		*scheme = 3;
		break;
	}

	switch (*scheme) {
	case 1:
		linearpred(pointsin, X, coeff, nbrs, remove, inter, nn, weights,
				docoeff);
		break;

	case 2:
		quadpred(pointsin, X, coeff, nbrs, remove, inter, nn, weights, docoeff);
		break;

	case 3:
		cubicpred(pointsin, X, coeff, nbrs, remove, inter, nn, weights,
				docoeff);
		break;
	}

}

/* ********** */

void amatdual(int *steps, int *po, int *re, int *nbrs, double *weights,
	double *alpha, int *lpo, int *lre, int *nn, double *adual)
{
	int n = *lpo + *lre, lnp = n - *steps + 1, l = lnp - 1;
	int i, j, k, m, r1, r2, ind, ind2, dummy;
	int *newpoints = calloc(l, sizeof(int)), *order = calloc(*nn, sizeof(int)),
			*revre = calloc(*lre, sizeof(int));

	double *malpha = calloc(*nn, sizeof(double)), *lastcol = calloc(l,
			sizeof(double)), *gdual = calloc(lnp, sizeof(double));
	double *hdual = calloc(l * lnp, sizeof(double)), *obyo = calloc(*nn * *nn,
			sizeof(double));
	double *tmp = calloc(l, sizeof(double)), *hdualtmp = calloc(l * l,
			sizeof(double));

	/*void mmult();
	void mycbind();
	void mydiag();
	void mymatchi();
	void myrbind();
	void myrevi();*/

	myrevi(re, lre, revre);

	for (i = 0; i < *lpo; i++)
		*(newpoints + i) = *po++;
	for (j = 0; j < (*lre - *steps); j++) /* *steps-1 is different from R code (Amatdual3) [see line 7 & 10] */
		*(newpoints + *lpo + j) = *(revre + j);

	for (k = 0; k < *nn; k++) {
		*(malpha + k) = -*(alpha + k);
	}

	free(revre);

	dummy = 1;

	mymatchi(nbrs, newpoints, nn, &l, order); /* finds out where nbrs come in hdual */

	free(newpoints);

	for (m = 0; m < *nn; m++) {
		*(lastcol + *(order + m) - 1) = *(alpha + m);
		*(gdual + *(order + m) - 1) = -*(weights + m);
	}

	mydiag(tmp, &l, &dummy, hdualtmp); /* non-specific diagonal in tmp */

	free(tmp);

	mycbind(hdualtmp, lastcol, &l, &l, &dummy, hdual);

	free(hdualtmp);
	free(lastcol);

	mmult(malpha, weights, nn, &dummy, nn, obyo);

	free(malpha);

	for (r1 = 0; r1 < *nn; r1++) {
		for (r2 = 0; r2 < *nn; r2++) {
			ind = ((*(order + r1) - 1) * lnp) + *(order + r2) - 1;
			ind2 = (*nn * r1) + r2;
			*(hdual + ind) = *(obyo + ind2);
		}
		*(hdual + (*(order + r1) - 1) * (lnp + 1)) += 1;
	}
	free(obyo);
	free(order);

	*(gdual + l) = 1;
	myrbind(hdual, gdual, &l, &lnp, &dummy, adual);

	free(gdual);
	free(hdual);
}

/* ********** */
/* does (pointwise) vector multiplication and i.p. */
void atimesb(double *a, double *b, int *n, double *prod, double *s)
{
	int i;
	*s = 0;

	for (i = 0; i < *n; i++) {
		*(prod + i) = *(a + i) * *(b + i);
		*s += *(prod + i);
	}

}

/* ********** */
/* fills the top left corner of Aaug with A.  If Aaug is initialized as zeros, this will augment A with a column and row of zeros.  */
/* input for A and Aaug are by row (in vector form). */
void aug(double *A, int *ra, int *ca, double *Aaug)
{
	int i = 1, j = *ca + 1;
	double *zeroc = calloc(*ra, sizeof(double)), *zeror = calloc(*ca + 1,
			sizeof(double));
	double *tmp = calloc(*ra * (*ca + 1), sizeof(double));

	/*void myrbind();
	void mycbind();*/

	mycbind(A, zeroc, ra, ca, &i, tmp);
	free(zeroc);
	myrbind(tmp, zeror, ra, &j, &i, Aaug);

	free(tmp);
	free(zeror);

}

/* ********** */

void cubicpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	int *inter, int *nn, double *weights, int *docoeff)
{
	int i, one = 1, two = 2, three = 3, nc = 3, dim;
	double *xn, *xr;
	double *tmpxn, *txn;
	double *tmp, *tmp2, *mm, *cn;
	double *tmpxr, *bhat;
	double *dummy1, *dummy2, *dummy3;
	double pred = 0;

	/*void mmult();
	void mycbind();
	//void mycpyd();
	void myt();
	void rmatsolve();*/

	tmpxr = calloc(1, sizeof(double));
	*tmpxr = *(X + *remove - 1);

	dummy1 = calloc(*nn, sizeof(double));
	dummy2 = calloc(*nn, sizeof(double));
	dummy3 = calloc(*nn, sizeof(double));
	cn = calloc(*nn, sizeof(double));

	for (i = 0; i < *nn; i++) {
		*(dummy1 + i) = *(X + *(nbrs + i) - 1);
		*(dummy2 + i) = pow(*(dummy1 + i), 2);
		*(dummy3 + i) = pow(*(dummy1 + i), 3);
		*(cn + i) = *(coeff + *(nbrs + i) - 1);
	}

	txn = calloc(MULTIPLY2(*nn), sizeof(double));
	tmpxn = calloc(MULTIPLY3(*nn), sizeof(double));
	mycbind(dummy1, dummy2, nn, &one, &one, txn);
	mycbind(txn, dummy3, nn, &two, &one, tmpxn);

	free(txn);
	free(dummy1);
	free(dummy2);
	free(dummy3);

	if (*inter == 1) {
		tmp = calloc(*nn, sizeof(double));
		for (i = 0; i < *nn; i++) {
			*(tmp + i) = 1.0;
		}
		xr = calloc(4, sizeof(double));
		xn = calloc(4 * *nn, sizeof(double));
		mycbind(tmp, tmpxn, nn, &one, &three, xn);
		*xr = 1.0;
		*(xr + 1) = *tmpxr;
		*(xr + 2) = pow(tmpxr[0], 2);
		*(xr + 3) = pow(tmpxr[0], 3);
		nc = 4;
		free(tmp);
	} else {
		xr = calloc(3, sizeof(double));
		xn = calloc(3 * *nn, sizeof(double));
		*xr = *tmpxr;
		*(xr + 1) = pow(tmpxr[0], 2);
		*(xr + 2) = pow(tmpxr[0], 3);
		dim = 3 * *nn;
		mycpyd(tmpxn, &dim, xn);
	}

	free(tmpxr);
	free(tmpxn);

	dim = nc * *nn;
	tmpxr = calloc(nc, sizeof(double));
	tmpxn = calloc(nc * *nn, sizeof(double));
	mycpyd(xr, &nc, tmpxr);
	mycpyd(xn, &dim, tmpxn);

	if (*nn == 3) {
		free(xn);
		xn = calloc((nc - 1) * *nn, sizeof(double));

		for (i = 0; i < (nc - 1); i++) {
			*(xn + i) = *(tmpxn + i);
			*(xn + nc + i - 1) = *(tmpxn + nc + i);
			*(xn + (2 * nc) + i - 2) = *(tmpxn + (2 * nc) + i);
		}
		nc -= 1;
		free(xr);
		xr = calloc(nc, sizeof(double));
		mycpyd(tmpxr, &nc, xr);
	}

	if (*nn == 2) {
		free(xn);
		xn = calloc((nc - 2) * *nn, sizeof(double));
		for (i = 0; i < (nc - 2); i++) {
			*(xn + i) = *(tmpxn + i);
			*(xn + nc + i - 2) = *(tmpxn + nc + i);
		}
		nc -= 2;
		free(xr);
		xr = calloc(nc, sizeof(double));
		mycpyd(tmpxr, &nc, xr);
	}

	free(tmpxn);
	free(tmpxr);

	if (*nn >= 2) {
		dim = nc * *nn;
		tmp = calloc(nc * nc, sizeof(double));
		txn = calloc(dim, sizeof(double));
		myt(xn, nn, &nc, txn);
		mmult(txn, xn, &nc, nn, &nc, tmp);
		tmp2 = calloc(nc * nc, sizeof(double));
		rmatsolve(tmp, &nc, tmp2);

		free(xn);
		free(tmp);

		mm = calloc(dim, sizeof(double));
		bhat = calloc(*nn, sizeof(double));
		mmult(tmp2, txn, &nc, &nc, nn, mm);
		mmult(mm, cn, &nc, nn, &one, bhat);

		free(txn);
		free(tmp2);
		free(cn);

		mmult(xr, bhat, &one, &nc, &one, &pred);
		mmult(xr, mm, &one, &nc, nn, weights);

		free(xr);
		free(mm);
		free(bhat);
	} else {
		*weights = 1;
		pred = *(coeff + *nbrs - 1);

		/*MAN: 2/4/09

		 xn,xr,cn still exist (even if not used) and need to be free'd for the 1 neighbour case: */

		free(xr);
		free(xn);
		free(cn);

	}

	if (*docoeff == 1) {
		*(coeff + *remove - 1) -= pred;
	}

}

/* ********** */

void fwtnp(double *input, double *f, int *nkeep, int *intercept,
	int *initboundhandl, int *neighbours, int *closest,int *LocalPred, int *n,
	double *coeff, double *lengthsremove,	double *lengths, double *lca,
	int *pointsin, int *nc, int *doW, double *W, int *varonly, double *v)
{

	int i, j, k = 0, N = *n, nn, scheme, r, remove, nr = 0, max, dim, dim1,
			dim2, dim2sq, dummy = *n - *nkeep, nnmax = MULTIPLY2(*neighbours), one = 1,
			ex = 0;
	int *po;
	int *nbrs2;
	int *index2;
	int *nbrs;
	int *index;

	double min;
	double *X = malloc(*n * sizeof(double));
	double *I = malloc((*n + 1) * sizeof(double));
	double *sX = malloc(*n * sizeof(double));

	double *weights2;
	double *len2;
	double *newline;
	double *tmplca;
	double *alpha, *weights;

	double *Wnew = 0, *Wtmp = 0;

	/*void adaptneigh();
	void adaptpred();
	void cubicpred();
	void delrow();
	void getnbrs();
	void getridd();
	void getridi();
	void intervals();
	void linearpred();
	void makelcaline();
	//void mycpyd();
	//void mycpyi();
	void mysortd();
	void mymind();
	void pointsupdate();
	void pts();
	void quadpred();
	void updatelca();*/

	mycpyd(input, n, X);
	intervals(X, initboundhandl, n, I);
	for (i = 0; i < *n; i++) {
		*(lengths + i) = *(I + i + 1) - *(I + i);
	}
	free(I);

	mysortd(X, n, sX, pointsin, &one);
	mycpyd(f, n, coeff);
	free(sX);

	/* (doW,varonly) should be (0,0),(1,0),or(0,1) */

	ex = *doW + *varonly;

	if (ex == 1) {
		Wnew = calloc(*n * *n, sizeof(double));
		for (i = 0; i < *n; i++) {
			*(Wnew + (i * *n) + i) = 1;
		}
	}
	/*	if(*varonly==1){
	 *(v+i)=1;
	 }*/

	if (*nkeep != *n) {
		for (j = 1; j <= dummy; j++) {
			mymind(lengths, &N, &min, &remove);
			remove = *(pointsin + remove - 1);

			nbrs = calloc(nnmax, sizeof(int)); /* set up as maximal */
			index = calloc(nnmax, sizeof(int)); /* ... */

			if (*LocalPred == 5) {

				nbrs2 = calloc(nnmax, sizeof(int));
				index2 = calloc(nnmax, sizeof(int));

				/*            mycpyi(nbrs,&nnmax,nbrs2);
				 mycpyi(index,&nnmax,index2); 	wierd?

				 mycpyi(nbrs,&nn,nbrs2);
				 mycpyi(index,&nn,index2);*/

				weights2 = calloc(nnmax, sizeof(double));

			} else { /* known nn given by getnbrs */
				getnbrs(X, &remove, pointsin, &N, neighbours, closest, nbrs,
						index, &nn);

				nbrs2 = calloc(nn, sizeof(int));
				index2 = calloc(nn, sizeof(int));

				mycpyi(nbrs, &nn, nbrs2); /*   fill to proper size */
				mycpyi(index, &nn, index2); /*   ... */

				weights2 = calloc(nn, sizeof(double));

			}
			free(nbrs);
			free(index);

			switch (*LocalPred) {
			case 1:
				scheme = 1;
				linearpred(pointsin, X, coeff, nbrs2, &remove, intercept, &nn,
						weights2, &one);
				break;

			case 2:
				scheme = 2;
				quadpred(pointsin, X, coeff, nbrs2, &remove, intercept, &nn,
						weights2, &one);
				break;

			case 3:
				scheme = 3;
				cubicpred(pointsin, X, coeff, nbrs2, &remove, intercept, &nn,
						weights2, &one);
				break;

			case 4:
				scheme = 1;
				adaptpred(pointsin, X, coeff, nbrs2, &remove, intercept, &nn,
						weights2, &scheme, &one);
				break;

			case 5:
				scheme = 1;
				adaptneigh(pointsin, X, coeff, nbrs2, &remove, intercept, &nn,
						weights2, &scheme, closest, index2, neighbours, &N,
						&one);
				break;
			}

			nbrs = calloc(nn, sizeof(int)); /* nn should be known in all cases now */
			index = calloc(nn, sizeof(int));
			weights = calloc(nn, sizeof(double));
			alpha = calloc(nn, sizeof(double));

			mycpyi(nbrs2, &nn, nbrs);
			mycpyi(index2, &nn, index);
			mycpyd(weights2, &nn, weights);

			free(nbrs2);
			free(index2);
			free(weights2);

			pointsupdate(X, coeff, &nn, index, &remove, pointsin, weights,
					lengths, &N, alpha, &r);

			*(lengthsremove + j - 1) = *(lengths + r - 1);

			newline = calloc((3 * nn + 5), sizeof(double));
			makelcaline(&remove, &nn, nbrs, alpha, weights, &scheme, intercept,
					closest, newline);
			max = (*nc >= (3 * nn + 5)) ? *nc : (3 * nn + 5);
			tmplca = calloc(max * j, sizeof(double));
			updatelca(lca, &nr, nc, newline, tmplca);

			if (ex == 1) {
				if (*varonly == 1) {
					for (i = 0; i < *n; i++) {
						for (k = 0; k < nn; k++) {
							*(Wnew + (r - 1) * *n + i) -= *(weights + k)
									* *(Wnew + (*(index + k) - 1) * *n + i);
						}
					}
					for (i = 0; i < *n; i++) {
						for (k = 0; k < nn; k++) {
							*(Wnew + (*(index + k) - 1) * *n + i) +=
									*(alpha + k) * *(Wnew + (r - 1) * *n + i);
						}
						*(v + remove - 1) += pow(*(Wnew + (r - 1) * *n + i), 2);
					}

					dim = *n - j + 1;
					dim1 = dim - 1;
					dim2 = dim * *n;
					Wtmp = calloc(dim2, sizeof(double));
					mycpyd(Wnew, &dim2, Wtmp);
					free(Wnew);
					Wnew = calloc(dim1 * *n, sizeof(double));
					delrow(Wtmp, &dim, n, &r, Wnew);
					free(Wtmp);
				} else {
					for (i = 0; i < *n; i++) {
						for (k = 0; k < nn; k++) {
							*(Wnew + (remove - 1) * *n + i) -= *(weights + k)
									* *(Wnew + (*(nbrs + k) - 1) * *n + i);
						}
					}
					for (i = 0; i < *n; i++) {
						for (k = 0; k < nn; k++) {
							*(Wnew + (*(nbrs + k) - 1) * *n + i) += *(alpha + k)
									* *(Wnew + (remove - 1) * *n + i);
						}
					}
				}
			}

			free(nbrs);
			free(alpha);
			free(weights);
			free(newline);
			free(index);

			len2 = calloc(N, sizeof(double));
			po = calloc(N, sizeof(int));

			mycpyd(lengths, &N, len2);
			mycpyi(pointsin, &N, po);

			getridd(len2, &N, &r, lengths);
			getridi(po, &N, &r, pointsin);
			free(len2);
			free(po);

			dim = nr * *nc;
			mycpyd(tmplca, &dim, lca);
			free(tmplca);
			N -= 1;
		} /* j */
	} /* if */

	if (ex == 1) {
		if (*varonly == 1) {
			for (i = 0; i < N; i++) {
				for (k = 0; k < *n; k++) {
					*(v + *(pointsin + i) - 1) += pow(*(Wnew + (i * *n) + k),
							2);
				}
			}
			dim1 = *nkeep * *n;
			mycpyd(Wnew, &dim1, W);
		} else {
			dim2sq = pow(*n, 2);
			mycpyd(Wnew, &dim2sq, W);
		}

		free(Wnew);
	}
	/* }        */

	free(X);
}

/* ********** */

void getnbrs2(double *X, int *remove, int *pointsin, int *lpo, int *neigh,
	int *closest, int *nbrs, int *index, int *nn)
{
	int d1, r, lr = (*lpo - 2), yn, i, j, ln = 0, rn = 0, lt, numleft = 0;
	int *range;
	int *checkr, *checkl;
	int *temp1, *q;
	int *temp2;
	int *tempq, *tempnbrs;
	int *tempindex, *B, *dummy, one = 1;
	double *sd, sl = 0, sr = 0, *distances;

	/*void mywhichi();
	void getridi();
	void mysortd();
	void mysorti();
	//void mycpyi();
	//void mycpyd();*/

	mywhichi(pointsin, lpo, remove, &r);

	range = calloc(lr, sizeof(int));
	for (j = 0; j < (*lpo - 2); j++) {
		*(range + j) = j + 2;
	}

	mywhichi(range, &lr, &r, &yn);
	/* yn= (yn==(lr+1)) ? 0 : 1;*/

	free(range);

	if (r == 1) {
		tempindex = calloc(1, sizeof(int));
		*tempindex = r + 1;
		tempnbrs = calloc(1, sizeof(int));
		*tempnbrs = *(pointsin + r);
		ln = 0;
		rn = 1;
	} else if (r == *lpo) {
		tempindex = calloc(1, sizeof(int));
		*tempindex = r - 1;
		tempnbrs = calloc(1, sizeof(int));
		*tempnbrs = *(pointsin + *tempindex - 1);
		ln = 1;
		rn = 0;
	}

	/*if(yn==1){*/
	else {
		checkr = calloc(*neigh, sizeof(int));
		checkl = calloc(*neigh, sizeof(int));
		for (i = 0; i < *neigh; i++) {
			*(checkl + i) = ((r - i - 1) < 1) ? 1 : 0; /* can avoid using checkr and checkl */
			*(checkr + i) = ((r + i + 1) > *lpo) ? 1 : 0;
			sl += *(checkl + i);
			sr += *(checkr + i);
		}
		ln = *neigh - (int) sl;
		rn = *neigh - (int) sr;
		free(checkl);
		free(checkr);
		lt = rn + ln;
		tempnbrs = calloc(lt, sizeof(int));
		tempindex = calloc(lt, sizeof(int)); /* gets valid neighbours according to initial request and pointsin */
		for (i = 0; i < ln; i++) {
			*(tempnbrs + i) = *(pointsin + r - ln + i - 1);
			*(tempindex + i) = r - ln + i;
		}
		for (i = 0; i < rn; i++) {
			*(tempnbrs + ln + i) = *(pointsin + r + i);
			*(tempindex + ln + i) = r + i + 1;
		}
	}

	if (*closest == 1) {
		distances = calloc(lt, sizeof(double));
		for (i = 0; i < lt; i++) {
			*(distances + i) = fabs(
					*(X + *remove - 1) - *(X + *(tempnbrs + i) - 1));
		}

		d1 = (*neigh <= (*lpo - 1)) ? *neigh : (*lpo - 1);
		sd = calloc(lt, sizeof(double));
		tempq = calloc(lt, sizeof(int));
		mysortd(distances, &lt, sd, tempq, &one);
		free(distances);
		free(sd);
		temp1 = calloc(lt, sizeof(int));
		temp2 = calloc(lt, sizeof(int));
		mycpyi(tempnbrs, &lt, temp1);
		mycpyi(tempindex, &lt, temp2);
		free(tempnbrs);
		free(tempindex);
		tempindex = calloc(d1, sizeof(int)); /* ******************* */
		tempnbrs = calloc(d1, sizeof(int));
		B = calloc(d1, sizeof(int));
		q = calloc(d1, sizeof(int));

		for (i = 0; i < d1; i++) {
			*(q + i) = *(tempq + i);
			*(tempnbrs + i) = *(temp1 + *(q + i) - 1);
			*(tempindex + i) = *(temp2 + *(q + i) - 1);
			*(B + i) = (*(tempindex + *(q + i) - 1) < r) ? 1 : 0; /* ****************** */
			numleft += *(B + i);
		}
		free(tempq);
		free(temp1);
		free(temp2);
		free(B);
		free(q);
		ln = numleft;
		rn = d1 - ln;
	}

	lt = ln + rn;
	dummy = calloc(lt, sizeof(int));
	temp2 = calloc(lt, sizeof(int));
	mycpyi(tempindex, &lt, temp2);
	mysorti(temp2, &lt, index, dummy, &one);
	free(dummy);
	free(temp2);
	for (i = 0; i < lt; i++) {
		*(nbrs + i) = *(pointsin + *(index + i) - 1);
	}
	*nn = lt;

}

/* ********** */

void getridd(double *a, int *la, int *pos, double *b) {
	int i;
	for (i = 0; i < (*pos - 1); i++) {
		*b++ = *a++;
	}
	a++;
	for (i = *pos; i < *la; i++) {
		*b++ = *a++;
	}

}

/* ********** */

void getridi(int *a, int *la, int *pos, int *b) {
	int i;
	for (i = 0; i < (*pos - 1); i++) {
		*b++ = *a++;
	}
	a++;
	for (i = *pos; i < *la; i++) {
		*b++ = *a++;
	}

}

/* ********** */

void intervals(double *X, int *initboundhandl, int *n, double *inter)
{
	double *sX = calloc(*n, sizeof(double));
	int i, one = 1;
	int l = *n, *order = calloc(*n, sizeof(int));

	//void mysortd();

	mysortd(X, n, sX, order, &one);

	free(order);

	switch (*initboundhandl) {
	case 0: /* reflect */
		*inter = (3 * (*sX) - *(sX + 1)) / 2;
		*(inter + l) = (3 * (*(sX + l - 1)) - *(sX + l - 2)) / 2;
		break;
	case 1: /* stop */
		*inter = *sX;
		*(inter + l) = *(sX + l - 1);
		break;
	} /* end switch */

	for (i = 0; i < (*n - 1); i++) {
		*(inter + i + 1) = (*(sX + i) + *(sX + i + 1)) / 2; /* can possibly use *inter++ somehow here? */
	}
	free(sX);
}

/* ********** */

void linearpred_nn_not_one(double *X, double *coeff, int *nbrs, int *remove,
	int *inter, int *nn, double *weights, double *tmp, double *pred){
	int one = 1, nc = 1, i;
	double *cn = malloc(*nn * sizeof(double)); //changed to malloc
	double *tmpxn = malloc(*nn * sizeof(double));  //changed to malloc

	for (i = 0; i < *nn; i++) {
		*(tmpxn + i) = *(X + *(nbrs + i) - 1);
		*(cn + i) = *(coeff + *(nbrs + i) - 1);
	}

	double tmpxr = *(X + *remove - 1);
	double *xn, *xr;
	if (*inter == 1) {
		tmp = malloc(*nn * sizeof(double));  //changed to malloc
		for (i = 0; i < *nn; i++) {
			*(tmp + i) = 1.0;
		}
		xn = malloc(MULTIPLY2(*nn) * sizeof(double)); //changed to malloc
		xr = malloc(2 * sizeof(double));  //changed to malloc
		mycbind(tmp, tmpxn, nn, &one, &one, xn);
		*xr = 1.0;
		*(xr + 1) = tmpxr;
		nc = 2;
		free(tmp);
	} else {
		xn = malloc(*nn * sizeof(double)); //changed to malloc
		xr = malloc(1 * sizeof(double)); //changed to malloc
		mycpyd(tmpxn, nn, xn);
		*xr = tmpxr;
	}
	free(tmpxn);

	int dim = nc * *nn;
	tmp = malloc(nc * nc * sizeof(double)); //changed to malloc
	double *txn = malloc(dim * sizeof(double)); //changed to malloc
	myt(xn, nn, &nc, txn);
	mmult(txn, xn, &nc, nn, &nc, tmp);
	double *tmp2 = malloc(nc * nc * sizeof(double)); //changed to malloc
	rmatsolve(tmp, &nc, tmp2);

	free(xn);
	free(tmp);


	double *mm = malloc(dim * sizeof(double)); //changed to malloc
	double *bhat = malloc(*nn * sizeof(double)); //changed to malloc
	mmult(tmp2, txn, &nc, &nc, nn, mm);
	mmult(mm, cn, &nc, nn, &one, bhat);

	free(txn);
	free(tmp2);
	free(cn);

	mmult(xr, bhat, &one, &nc, &one, pred);
	mmult(xr, mm, &one, &nc, nn, weights);
	free(xr);
	free(bhat);
	free(mm);
}

void linearpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	int *inter, int *nn, double *weights, int *docoeff)
{
	double pred, *tmp = 0;

	if (*nn == 1) {
		*weights = 1;
		pred = *(coeff + *nbrs - 1);
	} else {
		linearpred_nn_not_one(X, coeff, nbrs, remove, inter, nn, weights,
			tmp, &pred);
	}
	if (*docoeff == 1) {
		coeff[*remove - 1] -= pred;
		//*(coeff + *remove - 1) -= pred;
	}

}

/* ********** */
/* hopefully, this will coerce all necessary arguments to double and concatenate them */
void makelcaline(int *remove, int *nn, int *nbrs, double *alpha,
	double *weights, int *scheme, int *inter, int *closest,	double *newline)
{
	int i;

	*newline = (double) *remove;
	*(newline + 1) = (double) *nn;

	for (i = 0; i < *nn; i++) {
		*(newline + 2 + i) = (double) *(nbrs + i);
		*(newline + 2 + *nn + i) = *(alpha + i);
		*(newline + 2 + (2 * *nn) + i) = *(weights + i);
	}

	*(newline + 2 + MULTIPLY3(*nn)) = (double) *scheme;
	*(newline + 3 + MULTIPLY3(*nn)) = (double) *inter;
	*(newline + 4 + MULTIPLY3(*nn)) = (double) *closest;
}

/* ********** */
void mmult(double *x, double *y, int *nrx, int *ncx, int *ncy, double *ans)
{
	char *transa = "T", *transb = "T";
	double one = 1.0, zero = 0.0, *tmp = malloc(*nrx * *ncy * sizeof(double));

	//void myt();

	F77_CALL(dgemm)(transa, transb, nrx, ncy, ncx, &one, x, ncx, y, ncy, &zero,
			tmp, nrx);

	myt(tmp, ncy, nrx, ans); // output from call is bycol, i.e. need to swap dimensions

	free(tmp);
}

/* ********** */

void mycbind(double *a, double *b, int *ra, int *ca, int *cb, double *c)
{

	int i, j, k;

	for (i = 0; i < *ra; i++) {
		for (j = 0; j < *ca; j++) {
			*c++ = *a++;
		}
		for (k = 0; k < *cb; k++) {
			*c++ = *b++;
		}
	}

}

/* ********** */
//copy a in to b
void mycpyd(double *a, int *len, double *b)
{
	memcpy(b, a, *len * sizeof(double));
	/*int i;

	for (i = 0; i < *len; i++) {
		*(b + i) = *(a + i);
	}*/

}

/* **********
	a in to b
 */

void mycpyi(int *a, int *len, int *b)
{
	//void * memcpy ( void * destination, const void * source, size_t num );
	memcpy(b, a, *len * sizeof(int));
	/*int i = 0;
	for (i; i < *len; i++) {
		*(b + i) = *(a + i);
	}*/
}

/* ********** */

void mydiag(double *d, int *n, int *ones, double *m)
{

	int j, k;

	if (*ones == 1) { /* sets up the "default" */
		for (j = 0; j < *n; j++)
			*(d + j) = 1;
	}

	for (k = 0; k < (*n * *n); k++)
		*(m + k) = 0;

	for (k = 0; k < *n; k++)
		*(m + (k * *n) + k) = *(d + k);

}

/* ********** */

void mymatchd(double *numa, double *numb, int *lnuma, int *lnumb, int *pos)
{
	int i;

	void mywhichd();

	for (i = 0; i < *lnuma; i++) {
		mywhichd(numb, lnumb, numa, pos);
		numa++;
		pos++;
	}

}

/* ********** */

void mymatchi(int *numa, int *numb, int *lnuma, int *lnumb, int *pos)
{
	int i;

	void mywhichi();

	for (i = 0; i < *lnuma; i++) {
		mywhichi(numb, lnumb, numa, pos);
		numa++;
		pos++;
	}

}

/* ********** */

void mymaxd(double *num, int *lnum, double *max, int *pos)
{
	int i, start = 0;
	*max = *num;
	*pos = start;

	for (i = 0; i < (*lnum - 1); i++) {
		*pos = (*max >= *(num + 1)) ? *pos : (i + 1);
		*max = (*max >= *(num + 1)) ? *max : *(num + 1);
		num++;
	}
	*pos += 1;
}

/* ********** */

void mymaxi(int *num, int *lnum, int *max, int *pos)
{
	int i, start = 0;
	*max = *num;
	*pos = start;

	for (i = 0; i < (*lnum - 1); i++) {
		*pos = (*max >= *(num + 1)) ? *pos : (i + 1);
		*max = (*max >= *(num + 1)) ? *max : *(num + 1);
		num++;
	}
	*pos += 1;
}

/* ********** */

void mymind(double *num, int *lnum, double *min, int *pos)
{
	int i, start = 0;
	*min = *num;
	*pos = start;

	for (i = 0; i < (*lnum - 1); i++) {
		*pos = (*min <= *(num + 1)) ? *pos : (i + 1);
		*min = (*min <= *(num + 1)) ? *min : *(num + 1);
		num++;
	}
	*pos += 1;
}

/* ********** */

void mymini(int *num, int *lnum, int *min, int *pos)
{
	int i, start = 0;
	*min = *num;
	*pos = start;

	for (i = 0; i < (*lnum - 1); i++) {
		*pos = (*min <= *(num + 1)) ? *pos : (i + 1);
		*min = (*min <= *(num + 1)) ? *min : *(num + 1);
		num++;
	}
	*pos += 1;
}

/* ********** */

void myrbind(double *a, double *b, int *ra, int *ca, int *rb, double *c)
{

	int i, j, k;

	for (i = 0; i < *ra; i++) {
		for (j = 0; j < *ca; j++) {
			*c++ = *a++;
		}
	}
	for (k = 0; k < *rb; k++) {
		for (j = 0; j < *ca; j++) {
			*c++ = *b++;
		}
	}

}

/* ********** */

void myrevd(double *dx, int *n, double *dy)
{

	int incx = -1, incy = 1;

	F77_NAME(dcopy)(n, dx, &incx, dy, &incy);

}

/* ********** */

void myrevi(int *a, int *la, int *b)
{
	int *tmp = a + (*la - 1); /* *tmp points to the last element of "a" */
	int i;

	for (i = 0; i < *la; i++) {
		*b++ = *tmp--;
	}

}

/* ********** */

void mysortd(double *a, int *la, double *sorted, int *order, int *inc)
{
	int i, *o = calloc(*la, sizeof(int));
	double *s = calloc(*la, sizeof(double));

	//void mycpyd();
	//void mycpyi();
	/*void myrevd();
	void myrevi();*/

	for (i = 0; i < *la; i++) {
		*(o + i) = (i + 1);
	}

	mycpyd(a, la, s);
	rsort_with_index(s, o, *la);

	if (*inc == 0) {
		myrevi(o, la, order);
		myrevd(s, la, sorted);
	} else {
		mycpyi(o, la, order);
		mycpyd(s, la, sorted);
	}
	free(o);
	free(s);
}

/* ********** */

void mysorti2(int *a, int *la, int *sorted, int *order, int *inc)
{
	int i, j, newpos, oldpos = *la + 10, curlen = *la, *curleft = calloc(*la,
			sizeof(int)), *oldleft = calloc(*la, sizeof(int));
	int min;
	int *current = calloc(*la, sizeof(int)), *old = calloc(*la, sizeof(int));
	int *o = calloc(*la, sizeof(int));

	int *s = calloc(*la, sizeof(int));

	/*void mymini();
	//void mycpyi();
	void getridi();
	void myrevi();*/

	for (j = 0; j < *la; j++) {
		*(curleft + j) = j + 1;
	}

	mycpyi(a, la, current);
	mycpyi(curleft, la, oldleft);

	for (i = 0; i < *la; i++) {
		mymini(current, &curlen, &min, &newpos);
		*(s + i) = min;
		oldpos = newpos - 1;
		*(o + i) = *(curleft + oldpos);
		free(old);
		old = calloc(curlen, sizeof(int));
		mycpyi(current, &curlen, old);
		free(oldleft);
		oldleft = calloc(curlen, sizeof(int));
		mycpyi(curleft, &curlen, oldleft);
		free(current);
		free(curleft);
		current = calloc((curlen - 1), sizeof(int));
		curleft = calloc((curlen - 1), sizeof(int));
		getridi(old, &curlen, &newpos, current);
		getridi(oldleft, &curlen, &newpos, curleft);
		curlen -= 1;
	}
	free(old);
	free(oldleft);
	free(current);
	free(curleft);
	if (*inc == 0) {
		myrevi(o, la, order);
		myrevi(s, la, sorted);
	} else {
		mycpyi(o, la, order);
		mycpyi(s, la, sorted);
	}
	free(o);
	free(s);
}

/* ********** */

void mysvd(double *a, int *n, double *rvalues, double *rvectors,
	int *decreasing)
{
	//char jobv[1], range[1], uplo[1];
	int il, iu, m, itmp, info = 0;
	double vl = 0.0, vu = 0.0, abstol = 0.0, tmp;

	char jobv[] = {'V'};
	char range[] = {'A'};
	char uplo[] = {'L'};

	int lwork = -1;
	int liwork = -1;
	int *isuppz = (int *)malloc(MULTIPLY2(*n) * sizeof(int));
	double *rv = (double *)malloc(*n * sizeof(double));
	double *rz = (double *)malloc(*n * *n * sizeof(double));

	F77_CALL(dsyevr)(jobv, range, uplo, n, a, n, &vl, &vu, &il, &iu, &abstol,
			&m, rv, rz, n, isuppz, &tmp, &lwork, &itmp, &liwork, &info);

	lwork = (int) tmp;
	liwork = itmp;

	double *work = (double *)malloc(lwork * sizeof(double));
	int *iwork = (int *)malloc(liwork * sizeof(int));

	F77_CALL(dsyevr)(jobv, range, uplo, n, a, n, &vl, &vu, &il, &iu, &abstol,
			&m, rv, rz, n, isuppz, work, &lwork, iwork, &liwork, &info);

	work = realloc(work, *n * *n * sizeof(double));
	myt(rz, n, n, work);

	//MAN: added 2/4/09 to avoid mem. leakage
	free(iwork);
	free(isuppz);
	free(rz);

	if (*decreasing == 1) {
		int j;
		for (int i = 0; i < *n; i++) {
			for (j = 0; j < *n; j++) {
				*(rvectors + j + (*n * i)) = *(work + (*n * i) + (*n - j - 1));
			}
		}
		myrevd(rv, n, rvalues);
	} else {
		m = *n * *n;
		mycpyd(rv, n, rvalues);
		mycpyd(work, &m, rvectors);
	}

	//MAN: added 2/4/09 to avoid mem. leakage
	free(work);
	free(rv);
}

/* ********** */

void myt(double *a, int *ra, int *ca, double *ta)
{
	int i, j, indA = 0;

	for (j = 0; j < *ca; j++) {
		for (i = 0; i < *ra; i++) {
			indA = j + (*ca * i);
			*ta++ = a[indA];
		}
	}

}

/* ********** */

void mywhichd(double *num, int *lnum, double *a, int *pos)
{
	/*int start = 0;
	*pos = start;*/
	*pos = 0;

	for (int i = 0; i < *lnum; i++) {
		if (*num == *a) {
			break;
		} else {
			num++;
			*pos += 1;
		}
	}
	*pos += 1;

}

/* ********** */

void mywhichi(int *num, int *lnum, int *a, int *pos)
{
	/*int start = 0;
	*pos = start;*/
	*pos = 0;

	for (int i = 0; i < *lnum; i++) {
		if (*num == *a) {
			break;
		} else {
			num++;
			*pos += 1;
		}
	}
	*pos += 1;

}

/* ********** */
void pointsupdate(double *X, double *coeff, int *nn, int *index, int *remove,
	 int *pointsin, double *wts, double *l, int *N, double *alpha, int *r)
{

	int q1, i, j;
	double s = 0;
	//void mywhichi();

	mywhichi(pointsin, N, remove, r);

	if (*r >= 2 && *r <= (*N - 1)) {
		for (i = 0; i < *nn; i++) {
			*(l + *(index + i) - 1) += *(l + *r - 1) * *(wts + i);
		} /* end for */
	} /* end if */
	else {
		if (*r == 1)
			*(l + 1) += *l;
		else
			*(l + *N - 2) += *(l + *N - 1);
	} /* end else */

	if (*nn == 1) {
		*alpha = (*(l + *r - 1) / *(l + *index - 1));
		q1 = *(pointsin + *index - 1);
		*(coeff + q1 - 1) += (*alpha * *(coeff + *remove - 1));
	} /* end if */
	else {
		for (i = 0; i < *nn; i++) {
			s += pow(*(l + *(index + i) - 1), 2);
		} /* end for */
		for (j = 0; j < *nn; j++) {
			*(alpha + j) = (*(l + *r - 1) * *(l + *(index + j) - 1)) / s;
			q1 = *(pointsin + *(index + j) - 1);
			*(coeff + q1 - 1) += *(alpha + j) * *(coeff + *remove - 1);
		} /* end for */
	} /* end else */

}

/* ********** */

void pts(double *input, double *start, int *n, double *X)
{
	double *ptr = &input[0];
	double *ptr2 = &X[0];
	X[0] = *start;

	for (int i = 1; i <= *n; i++){
		X[i] = *ptr2++ + *ptr++;
	}

}

/* ********** */

void quadpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	int *inter, int *nn, double *weights, int *docoeff)
{
	int i, one = 1, two = 2, nc = 2, dim;
	double *xn, *xr;
	double *tmpxn, *txn;
	double *tmp, *tmp2, *mm, *cn;
	double *tmpxr, *bhat, *dummy1, *dummy2;
	double pred = 0;

	/*void mmult();
	void mycbind();
	//void mycpyd();
	void myt();
	void rmatsolve();*/

	tmpxr = calloc(1, sizeof(double));
	*tmpxr = *(X + *remove - 1);

	dummy1 = calloc(*nn, sizeof(double));
	dummy2 = calloc(*nn, sizeof(double));
	cn = calloc(*nn, sizeof(double));

	for (i = 0; i < *nn; i++) {
		*(dummy1 + i) = *(X + *(nbrs + i) - 1);
		*(dummy2 + i) = pow(*(dummy1 + i), 2);
		*(cn + i) = *(coeff + *(nbrs + i) - 1);
	}

	tmpxn = calloc(MULTIPLY2(*nn), sizeof(double));

	mycbind(dummy1, dummy2, nn, &one, &one, tmpxn);

	free(dummy1);
	free(dummy2);

	if (*inter == 1) {
		tmp = calloc(*nn, sizeof(double));
		for (i = 0; i < *nn; i++) {
			*(tmp + i) = 1.0;
		}
		xr = calloc(3, sizeof(double));
		xn = calloc(3 * *nn, sizeof(double));
		mycbind(tmp, tmpxn, nn, &one, &two, xn);
		*xr = 1.0;
		*(xr + 1) = *tmpxr;
		*(xr + 2) = pow(tmpxr[0], 2);
		free(tmp);
		nc = 3;
	} else {
		xr = calloc(2, sizeof(double));
		xn = calloc(2 * *nn, sizeof(double));
		*xr = *tmpxr;
		*(xr + 1) = pow(tmpxr[0], 2);
		dim = MULTIPLY2(*nn);
		mycpyd(tmpxn, &dim, xn);
	}

	free(tmpxn);
	free(tmpxr);

	dim = nc * *nn;

	if (*nn == 2) {
		tmpxr = calloc(nc, sizeof(double));
		tmpxn = calloc(nc * *nn, sizeof(double));
		mycpyd(xr, &nc, tmpxr);
		mycpyd(xn, &dim, tmpxn);

		free(xn);
		xn = calloc((nc - 1) * *nn, sizeof(double));
		for (i = 0; i < (nc - 1); i++) {
			*(xn + i) = *(tmpxn + i);
			*(xn + nc + i - 1) = *(tmpxn + nc + i);
		}
		nc -= 1;
		free(xr);
		xr = calloc(nc, sizeof(double));
		mycpyd(tmpxr, &nc, xr);
		free(tmpxn);
		free(tmpxr);
	}

	if (*nn >= 2) {
		dim = nc * *nn;
		tmp = calloc(nc * nc, sizeof(double));
		txn = calloc(dim, sizeof(double));
		myt(xn, nn, &nc, txn);
		mmult(txn, xn, &nc, nn, &nc, tmp);
		tmp2 = calloc(nc * nc, sizeof(double));
		rmatsolve(tmp, &nc, tmp2);

		free(xn);
		free(tmp);

		mm = calloc(dim, sizeof(double));
		bhat = calloc(*nn, sizeof(double));
		mmult(tmp2, txn, &nc, &nc, nn, mm);
		mmult(mm, cn, &nc, nn, &one, bhat);

		free(txn);
		free(tmp2);
		free(cn);

		mmult(xr, bhat, &one, &nc, &one, &pred);
		mmult(xr, mm, &one, &nc, nn, weights);

		free(xr);
		free(mm);
		free(bhat);
	} else {
		*weights = 1;
		pred = *(coeff + *nbrs - 1);

		/*MAN: 2/4/09

		 xn,xr,cn still exist (even if not used) and need to be free'd for the 1 neighbour case: */

		free(xr);
		free(xn);
		free(cn);

	}
	if (*docoeff == 1) {
		*(coeff + *remove - 1) -= pred;
	}

}

/* ********** */

void rmatsolve(double *m, int *n, double *inv)
{
	double tmpd;

	if (*n == 1) {
		tmpd = *m;
		*inv = 1 / tmpd;
	}
	else {
		int zero = 0, one = 1;
		double *tmp, *d, *tev, *D, *rvalues, *rvectors;
		rvalues = malloc(*n * sizeof(double));
		rvectors = malloc(*n * *n * sizeof(double));
		d = malloc(*n * sizeof(double));

		mysvd(m, n, rvalues, rvectors, &one);
		for (int i = 0; i < *n; i++) {
			tmpd = *(rvalues + i);
			*(d + i) = 1 / tmpd;
		}
		free(rvalues);
		D = malloc(*n * *n * sizeof(double));
		tev = malloc(*n * *n * sizeof(double));
		tmp = malloc(*n * *n * sizeof(double));
		mydiag(d, n, &zero, D);
		free(d);
		myt(rvectors, n, n, tev);
		mmult(rvectors, D, n, n, n, tmp);
		free(D);
		free(rvectors);
		mmult(tmp, tev, n, n, n, inv);
		free(tmp);
		free(tev);
	}
}

/* ********** */
/* takes in the fwtnp lifting coefficient array, of size length(removelist)*/
/*by 3*max(n_r)+5 (zero filled,remove,nn,nbrs,alpha,weights,scheme,int,closest)*/
/*matno is length(rem), nc is ncol(lca) */
void transmatdual(double *lca, int *po, int *matno, int *lpo, int *nc,
	double *W, int *re)
{
	int i, j, steps, matsize, matdim, tolddim, nn, *nbrs;
	double *alpha, *weights, *A, *augment;

	/*void amatdual();
	void aug();
	void mmult();*/

	/*set up first matrices... */
	nn = (int) *(lca + (*nc * (*matno - 1)) + 1);

	nbrs = calloc(nn, sizeof(int));
	alpha = calloc(nn, sizeof(double));
	weights = calloc(nn, sizeof(double));

	for (i = 0; i < nn; i++) {
		*(nbrs + i) = (int) *(lca + (*nc * (*matno - 1)) + 2 + i);
		*(alpha + i) = *(lca + (*nc * (*matno - 1)) + nn + 2 + i);
		*(weights + i) = *(lca + (*nc * (*matno - 1)) + (2 * nn) + 2 + i);
	}

	amatdual(matno, po, re, nbrs, weights, alpha, lpo, matno, &nn, W);

	free(nbrs);
	free(alpha);
	free(weights);

	/*printf("done initial matrices\n");*/

	if (*matno > 1) {
		for (j = 2; j <= *matno; j++) {
			printf("j:%d\n", j);

			matdim = *lpo + j;
			tolddim = matdim - 1;
			matsize = matdim * matdim;
			steps = *matno - j + 1;

			augment = calloc(matsize, sizeof(double));
			aug(W, &tolddim, &tolddim, augment);
			*(augment + matsize - 1) = 1;

			nn = (int) *(lca + (*nc * (*matno - j)) + 1);
			nbrs = calloc(nn, sizeof(int));
			alpha = calloc(nn, sizeof(double));
			weights = calloc(nn, sizeof(double));

			for (i = 0; i < nn; i++) {
				*(nbrs + i) = (int) *(lca + (*nc * (*matno - j)) + 2 + i);
				*(alpha + i) = *(lca + (*nc * (*matno - j)) + nn + 2 + i);
				*(weights + i) =
						*(lca + (*nc * (*matno - j)) + (2 * nn) + 2 + i);
			}

			A = calloc(matsize, sizeof(double));
			amatdual(&steps, po, re, nbrs, weights, alpha, lpo, matno, &nn, A);
			free(nbrs);
			free(alpha);
			free(weights);
			mmult(augment, A, &matdim, &matdim, &matdim, W);
			free(augment);
			free(A);
		}
	}

}

/* ********** */

void undopointsupdate(double *coeff, int *nbrs, int *index, int *remove, int *r,
	int *N, double *gamweights, double *l, double *lr, double *alpha, int *nn)
{
	int j;
	double pred = 0, s = 0;

	if (*r >= 2 && *r <= *N) {
		for (j = 0; j < *nn; j++) {
			s += pow(*(l + *(index + j) - 1), 2);
		}
		for (j = 0; j < *nn; j++) {
			*(alpha + j) = (*(l + *(index + j) - 1) * *lr) / s;
			*(coeff + *(nbrs + j) - 1) -= *(alpha + j) * *(coeff + *remove - 1);
			*(l + *(index + j) - 1) -= *(gamweights + j) * *lr;
			pred += *(gamweights + j) * *(coeff + *(nbrs + j) - 1);
		}
	} else {
		*alpha = *lr / *(l + *index - 1);
		*(coeff + *nbrs - 1) -= (*alpha * *(coeff + *remove - 1));
		*(l + *index - 1) -= *lr;
		pred = *(coeff + *nbrs - 1);
	}

	*(coeff + *remove - 1) += pred;
}

/* ********** */
/* adds newline to lca, adding zeros where appropriate. */
/* NOTE: ncol(lca) is (3*nmax+5). */
/* initial nr is nrow(lca) */
/* initial nmax is ncol(lca) */
void updatelca(double *lca, int *nr, int *nc, double *newline, double *newlca)
{
	int one = 1, nn, ln, diff, lcadim, newnc;
	double *zeroadd, *tmplca, *tmpnl;

	/*void myrbind();
	void mycbind();
	//void mycpyd();*/

	nn = (int) *(newline + 1); /* new no. nbrs */
	ln = 3 * nn + 5; /* length(newline) */
	lcadim = *nc * *nr; /* present dim(lca) */

	if (lcadim == 0) {
		mycpyd(newline, &ln, newlca);
		*nr = 1;
		*nc = ln;
	} else {
		newnc = (*nc >= ln) ? *nc : ln;
		diff = ln - *nc;

		tmplca = calloc((*nr * newnc), sizeof(double)); /* extended lca */
		tmpnl = calloc(newnc, sizeof(double)); /* new newline */

		if (diff > 0) { /* ln is longer */
			zeroadd = calloc(diff * *nr, sizeof(double));
			mycbind(lca, zeroadd, nr, nc, &diff, tmplca);
			mycpyd(newline, &ln, tmpnl); /* ln is newnc */
			free(zeroadd);
		}
		if (diff < 0) { /* nc is longer */
			diff = -diff;
			zeroadd = calloc(diff, sizeof(double));
			mycbind(newline, zeroadd, &one, &ln, &diff, tmpnl);
			mycpyd(lca, &lcadim, tmplca); /* lca is already of dim (nr, newnc) */
			free(zeroadd);
		}
		if (diff == 0) { /* no difference, so do nothing, just add line later */
			mycpyd(lca, &lcadim, tmplca); /* new dim is old dim */
			mycpyd(newline, &ln, tmpnl);
		}

		*nc = newnc; /* update *nc to be new no. cols */

		myrbind(tmplca, tmpnl, nr, nc, &one, newlca); /* add newline */
		(*nr)++; /* update no. rows */

		free(tmplca);
		free(tmpnl);

	}

}

/* ********** */
/*gets schemehist from lca:

 1: LP
 2: QP
 3: CP

 */
void schfromlca(double *lca, int *nr, int *nc, int *sch)
{
	int i, nn;

	for (i = 0; i < *nr; i++) {
		nn = (int) *(lca + (i * *nc) + 1);
		*(sch + i) = (int) *(lca + (i * *nc) + (3 * nn + 2));
	}

}

/* ********** */
/*gets interhist from lca:

 0: no intercept
 1: intercept

 */
void interfromlca(double *lca, int *nr, int *nc, int *inter)
{
	int i, nn;

	for (i = 0; i < *nr; i++) {
		nn = (int) *(lca + (i * *nc) + 1);
		*(inter + i) = (int) *(lca + (i * *nc) + (3 * nn + 3));
	}

}

/* ********** */
/*gets clohist from lca:

 0: closest=FALSE
 1: closest=TRUE

 */
void clofromlca(double *lca, int *nr, int *nc, int *clo)
{
	int i, nn;

	for (i = 0; i < *nr; i++) {
		nn = (int) *(lca + (i * *nc) + 1);
		*(clo + i) = (int) *(lca + (i * *nc) + (3 * nn + 4));
	}

}

/* row no is given as R index (C index +1) */
void nbrsfromlca(double *lca, int *nc, int *rowno, int *nbrs)
{
	int i, nn;

	nn = *(lca + ((*rowno - 1) * *nc) + 1);

	for (i = 0; i < nn; i++) {
		*(nbrs + i) = (int) *(lca + ((*rowno - 1) * *nc) + 2 + i);
	}

}

/* row no is given as R index (C index +1) */
void afromlca(double *lca, int *nc, int *rowno, double *alpha)
{
	int i, nn;

	nn = *(lca + ((*rowno - 1) * *nc) + 1);

	for (i = 0; i < nn; i++) {
		*(alpha + i) = *(lca + ((*rowno - 1) * *nc) + 2 + nn + i);
	}

}

/* row no is given as R index (C index +1) */
void wfromlca(double *lca, int *nc, int *rowno, double *weights)
{
	int i, nn;

	nn = *(lca + ((*rowno - 1) * *nc) + 1);

	for (i = 0; i < nn; i++) {
		*(weights + i) = *(lca + ((*rowno - 1) * *nc) + 2 + (2 * nn) + i);
	}

}

void invtnp(double *X, double *coeff, double *lengths, double *lengthsremove,
	int *pointsin, double *lca, int *nadd, int *N, int *lr, int *nc, int *outpo,
	double *outlen)
{

	double lrem, *alpha, *weights, *dtmp;
	int n = *lr + *N, j, k, remove, nn, *nbrs, lcarow, *index, r, *itmp,
			one = 1, zero = 0;

	double *sX;
	int *o, *o2, Nplus, *q, *q1;
	int *schlist, *interlist;
	int scheme, inter;

	/*void schfromlca();
	void clofromlca();
	void interfromlca();
	//void mycpyd();
	//void mycpyi();
	void mymatchi();
	void mywhichi();
	void nbrsfromlca();
	void undopointsupdate();
	void mysortd();
	void mysorti();
	void linearpred();
	void quadpred();
	void cubicpred();*/

	mycpyi(pointsin, N, outpo); /* just in case they are not initialized as so */
	mycpyd(lengths, N, outlen); /* ... */

	sX = calloc(n, sizeof(double));
	o = calloc(n, sizeof(int));
	mysortd(X, &n, sX, o, &one);

	free(sX);

	schlist = calloc(*lr, sizeof(int));
	interlist = calloc(*lr, sizeof(int));

	schfromlca(lca, lr, nc, schlist);
	interfromlca(lca, lr, nc, interlist);

	if (*nadd > 0) {
		for (j = 1; j <= *nadd; j++) {
			remove = (int) *(lca + (*lr - j) * *nc);
			lrem = *(lengthsremove + (*lr - j));
			nn = (int) *(lca + (*lr - j) * *nc + 1);

			nbrs = calloc(nn, sizeof(int));
			index = calloc(nn, sizeof(int));

			lcarow = *lr - j + 1;
			nbrsfromlca(lca, nc, &lcarow, nbrs);
			mymatchi(nbrs, outpo, &nn, N, index);

			o2 = calloc(*N + 1, sizeof(int));
			q = calloc(*N + 1, sizeof(int));
			q1 = calloc(*N + 1, sizeof(int));
			Nplus = *N + 1;

			mycpyi(outpo, N, o2);
			*(o2 + *N) = remove;
			mymatchi(o2, o, &Nplus, &n, q1);
			mysorti(q1, &Nplus, o2, q, &one);
			mywhichi(q, &Nplus, &Nplus, &r); /* r should now be the (R) index of the newpoint [remove] into the new outpo */

			free(o2);
			free(q);
			free(q1);

			alpha = calloc(nn, sizeof(double));
			weights = calloc(nn, sizeof(double));

			scheme = *(schlist + *lr - j);
			inter = *(interlist + *lr - j);

			switch (scheme) { /* docoeff = 0 so that it doesn't affect coeff (just gets weights) */
			case 1:
				linearpred(pointsin, X, coeff, nbrs, &remove, &inter, &nn,
						weights, &zero);
				break;

			case 2:
				quadpred(pointsin, X, coeff, nbrs, &remove, &inter, &nn,
						weights, &zero);
				break;

			case 3:
				cubicpred(pointsin, X, coeff, nbrs, &remove, &inter, &nn,
						weights, &zero);
				break;
			}

			undopointsupdate(coeff, nbrs, index, &remove, &r, N, weights,
					outlen, &lrem, alpha, &nn);

			free(nbrs);
			free(index);
			free(alpha);
			free(weights);

			dtmp = calloc(*N, sizeof(double));
			itmp = calloc(*N, sizeof(int));
			mycpyd(outlen, N, dtmp);
			mycpyi(outpo, N, itmp);

			if (r == 1) {
				mycbind(&lrem, dtmp, &one, &one, N, outlen);
				*outpo = remove;
				for (k = 0; k < *N; k++) {
					*(outpo + k + 1) = *(itmp + k); /* mycbind is only for doubles */
				}
			} else if (r == (*N + 1)) {
				mycbind(dtmp, &lrem, &one, N, &one, outlen);
				mycpyi(itmp, N, outpo); /* mycbind is only for doubles */
				*(outpo + *N) = remove;
			} else {
				for (k = 0; k < (r - 1); k++) {
					*(outlen + k) = *(dtmp + k);
					*(outpo + k) = *(itmp + k);
				}
				*(outlen + r - 1) = lrem;
				*(outpo + r - 1) = remove;
				for (k = (r - 1); k < *N; k++) {
					*(outlen + k + 1) = *(dtmp + k);
					*(outpo + k + 1) = *(itmp + k);
				}
			}
			free(dtmp);
			free(itmp);
			*N += 1;
		}

	}

	free(o);
	free(schlist);
	free(interlist);
}

void findadds(int *rem, int *l, double *lca, int *nc, int *index, int *li,
	int *a)
{

	int *init2 = calloc(*li, sizeof(int));
	int *init = calloc(*li, sizeof(int));
	int *tmp2 = calloc(*li, sizeof(int));
	int *tmp, *nbrs, nn;
	int i, ii, j, jp, k, m, p;

	/*void mymatchi();
	void nbrsfromlca();
	void mywhichi();
	void mymaxi();*/
	//void mycpyi();

	tmp = calloc(*li, sizeof(int));

	mymatchi(index, rem, li, l, tmp);

	for (i = 0; i < *li; i++) {
		*(init + i) = *l - (*(tmp + i) - 1);
	}
	free(tmp);
	if (*l == 0) {
		for (i = 0; i < *li; i++) {
			*(init2 + i) = 1;
		}

	} else {
		for (i = 0; i < *li; i++) {
			tmp = calloc(*l, sizeof(int));
			m = 0;
			ii = *(index + i);
			for (j = 0; j < *l; j++) {
				jp = j + 1;
				nn = *(lca + (*nc * j) + 1);
				nbrs = calloc(nn, sizeof(int));
				nbrsfromlca(lca, nc, &jp, nbrs);
				mywhichi(nbrs, &nn, &ii, &k);
				*(tmp + j) = (k == (nn + 1)) ? 0 : jp;
				if (*(tmp + j) > 0) {
					m = jp;
					break;
				}
				free(nbrs);
			}
			*(init2 + i) = (m == 0) ? (*l + 1) : m;
			free(tmp);
		}
	}

	for (i = 0; i < *li; i++) {
		*(tmp2 + i) = *l - (*(init2 + i) - 1);
	}
	mycpyi(tmp2, li, init2);

	for (i = 0; i < *li; i++) {
		*(tmp2 + i) = (*(init + i) > *(init2 + i)) ? *(init + i) : *(init2 + i);
	}
	free(init);
	free(init2);

	mymaxi(tmp2, li, a, &p);

	free(tmp2);
}

void delrow(double *M, int *nr, int *nc, int *i, double *Mnew)
{

	int j, k, nrn = *nr - 1, dim = *nc * nrn;

	double *Mtmp = calloc(nrn * *nc, sizeof(double));

	for (j = 0; j < (*i - 1); j++) {
		for (k = 0; k < *nc; k++) {
			*(Mtmp + (j * *nc) + k) = *(M + (j * *nc) + k);
			/*		*Mtmp++=*(M+(j**nc)+k);*/
		}
	}
	for (j = *i; j < *nr; j++) {
		for (k = 0; k < *nc; k++) {
			*(Mtmp + ((j - 1) * *nc) + k) = *(M + (j * *nc) + k);
			/*		*Mtmp++=*(M+(j**nc)+k);*/
		}
	}

	mycpyd(Mtmp, &dim, Mnew);
	free(Mtmp);
}

void getnbrscase3(double *X, int *remove, int *pointsin, int *lpo, int *neigh,
	int *closest, int *nbrs, int *index, int *nn, int *r){
	/* find which exist out of neighbours on either side: */

	int i, left = 0, right = 0, tmp, n1;
	int *rr, dummy, one = 1;
	int *q, *tempindex;
	double *sd, *distances;

	for (i = 0; i < *neigh; i++) {
		tmp = ((*r - i - 1) > 0) ? 1 : 0;
		left += tmp;
		tmp = ((*r + i + 1) < (*lpo + 1)) ? 1 : 0;
		right += tmp;
	}
	tmp = left + right;
	rr = malloc(tmp * sizeof(int)); //changed to malloc
	for (i = 0; i < left; i++) {
		*(rr + i) = *r - left + i;
	}
	for (i = 0; i < right; i++) {
		*(rr + left + i) = *r + i + 1;
	}
	if (*closest == 1) {
		distances = malloc(tmp * sizeof(double)); //changed to malloc
		for (i = 0; i < tmp; i++) {
			dummy = *(pointsin + *(rr + i) - 1);
			*(distances + i) = fabs(*(X + dummy - 1) - *(X + *remove - 1));
		}
		sd = malloc(tmp * sizeof(double));  //changed to malloc
		q = malloc(tmp * sizeof(int));  //changed to malloc
		mysortd(distances, &tmp, sd, q, &one);
		free(distances);
		free(sd);
		n1 = (*neigh <= (*lpo - 1)) ? *neigh : (*lpo - 1);
		tempindex = calloc(n1, sizeof(int));
		for (i = 0; i < n1; i++) {
			*(tempindex + i) = *(rr + *(q + i) - 1);

		}
		free(q);
		q = malloc(n1 * sizeof(int));  //changed to malloc
		mysorti(tempindex, &n1, index, q, &one);
		free(q);
		free(tempindex);
		tmp = n1;
	}
	else { /* index should be original checked ones */
		for (i = 0; i < tmp; i++) {
			*(index + i) = *(rr + i);
		}
	}

	/*MAN: 2/4/09.  rr is used in case 3 whether closest or not, so should free outside if */
	free(rr);

	for (i = 0; i < tmp; i++) {
		*(nbrs + i) = *(pointsin + *(index + i) - 1);
	}
	*nn = tmp; /* tmp should exist (closest or not) */
}

void getnbrs(double *X, int *remove, int *pointsin, int *lpo, int *neigh,
	int *closest, int *nbrs, int *index, int *nn)
{

	int r, st;

	mywhichi(pointsin, lpo, remove, &r);
	if (r == 1) {
		st = 1;
	}
	if (r == *lpo) {
		st = 2;
	}
	if ((r > 1) & (r < *lpo)) { /* in "range" */
		st = 3;
	}

	switch (st) {
	case 1:
		*nbrs = *(pointsin + 1);
		*index = 2;
		*nn = 1;
		break;

	case 2:
		*nbrs = *(pointsin + *lpo - 2);
		*index = *lpo - 1;
		*nn = 1;
		break;

	case 3:
		getnbrscase3(X, remove, pointsin, lpo, neigh, closest, nbrs, index, nn, &r);
		break;
	} /* end switch  */
}

void mysorti(int *a, int *la, int *sorted, int *order, int *inc)
{

	int i;
	int *o = calloc(*la, sizeof(int)), *s2 = calloc(*la, sizeof(int));
	double *s = calloc(*la, sizeof(double));

	for (i = 0; i < *la; i++) {
		*(o + i) = i + 1;
		*(s + i) = (double) *(a + i);
	}

	/* this is a brute force, slightly inefficient way of
	 sorting integers with index since nothing else seems to work */

	//void mycpyi();
	//void myrevi();

	mycpyi(a, la, s2);

	rsort_with_index(s, o, *la); /* does the order bit... */

	R_isort(s2, *la); /* and now the sorting */

	if (*inc == 0) {
		myrevi(o, la, order);
		myrevi(s2, la, sorted);
	} else {
		mycpyi(o, la, order);
		mycpyi(s2, la, sorted);
	}
	free(o);
	free(s);
	free(s2);

}

/***************/

/*  "OTHER" VERSIONS OF SOME FUNCTIONS (INEFFICIENT)....... */

/*void mmult4(A,B,ra,ca,cb,C)

 does matrix multiplication.  A,B,C are inputted as vectors, with byrow=T

 double *A,*B,*C;
 int *ra,*ca,*cb;
 {
 double sum=0;
 int i,j,idxA, endA, idxB;

 for(i=0;i<*ra;i++){
 for(j=0;j<*cb;j++){
 idxA = i* *ca;
 endA = idxA + *ca;
 idxB = j;
 sum = 0;
 while (idxA<endA){
 sum += A[idxA] * B[idxB];
 idxA++;
 idxB += *cb;
 }
 *C++ = sum;
 }
 }

 }
 */
/*
 void myrevd(a,la,b)

 double *a,*b;
 int *la;
 {

 double *tmp=a+(*la-1);
 int i;

 for(i=0;i<*la;i++){
 *b++=*tmp--;
 }

 }

 void mysortd(a,la,sorted,order)

 double *a,*sorted;
 int *la, *order;
 {
 int i,j,newpos,oldpos=*la+10,curlen=*la,*curleft=malloc(*la*sizeof(int)),*oldleft=malloc(*la*sizeof(int));
 double min;
 double *current=malloc(*la*sizeof(double)),*old=malloc(*la*sizeof(double));

 void mymind();
 //void mycpyd();
 //void mycpyi();
 void getridd();
 void getridi();

 for(j=0;j<*la;j++){
 *(curleft+j)=j+1;
 }

 mycpyd(a,la,current);
 mycpyi(curleft,la,oldleft);

 for(i=0;i<*la;i++){
 mymind(current,&curlen,&min,&newpos);
 *(sorted+i)=min;
 oldpos=newpos-1;
 *(order+i)=*(curleft+oldpos);
 realloc(old,curlen*sizeof(double));
 mycpyd(current,&curlen,old);
 realloc(oldleft,curlen*sizeof(int));
 mycpyi(curleft,&curlen,oldleft);
 realloc(current,(curlen-1)*sizeof(double));
 realloc(curleft,(curlen-1)*sizeof(int));
 getridd(old,&curlen,&newpos,current);
 getridi(oldleft,&curlen,&newpos,curleft);
 curlen-=1;
 }

 }


 */

void mmult2(double *A, double *B, int *ra, int *ca, int *cb, double *C)
{
	double sum = 0;
	int i, j, k;

	for (i = 0; i < *ra; i++) {
		for (j = 0; j < *cb; j++) {
			sum = 0;
			for (k = 0; k < *ca; k++) {
				sum += *(A + (i * *ca) + k) * *(B + (k * *cb) + j);
			}
			*(C + (i * *cb) + j) = sum;
			/*printf("%lf",sum);*/
		}
	}

}

void mmult3(double *A, double *B, int *ra, int *ca, int *cb, double *C)
{
	double sum = 0;
	int i, j, k;

	int ra1 = *ra, ca1 = *ca, cb1 = *cb;

	for (i = 0; i < ra1; i++) {
		for (j = 0; j < cb1; j++) {
			sum = 0;
			for (k = 0; k < ca1; k++) {
				sum += A[(i * ca1) + k] * B[(k * cb1) + j];
			}
			C[(i * cb1) + j] = sum;
			/*printf("%lf",sum);*/
		}
	}

}
