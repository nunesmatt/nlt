#ifndef __ADLIFT
#define __ADLIFT

void adaptneigh(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	 int *inter, int *nn, double *weights, int *scheme, int *clo, int *index,
	 int *neighbours, int *N, int *docoeff);

void adaptpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	 int *inter, int *nn, double *weights, int *scheme, int *docoeff);

void amatdual(int *steps, int *po, int *re, int *nbrs, double *weights,
	double *alpha, int *lpo, int *lre, int *nn, double *adual);

/**
 * does (pointwise) vector multiplication and i.p.
 **/
void atimesb(double *a, double *b, int *n, double *prod, double *s);

/**
 * fills the top left corner of Aaug with A.  If Aaug is initialized as zeros,
 * this will augment A with a column and row of zeros. input for A and Aaug
 * are by row (in vector form).
 **/
void aug(double *A, int *ra, int *ca, double *Aaug);

void cubicpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	int *inter, int *nn, double *weights, int *docoeff);

void fwtnp(double *input, double *f, int *nkeep, int *intercept,
	int *initboundhandl, int *neighbours, int *closest,int *LocalPred, int *n,
	double *coeff, double *lengthsremove,	double *lengths, double *lca,
	int *pointsin, int *nc, int *doW, double *W, int *varonly, double *v);

void getnbrs2(double *X, int *remove, int *pointsin, int *lpo, int *neigh,
	int *closest, int *nbrs, int *index, int *nn);

void getridd(double *a, int *la, int *pos, double *b);
void getridi(int *a, int *la, int *pos, int *b);

void intervals(double *X, int *initboundhandl, int *n, double *inter);

void linearpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	int *inter, int *nn, double *weights, int *docoeff);

/**
 *hopefully, this will coerce all necessary arguments to double and
 *concatenate them
 **/
void makelcaline(int *remove, int *nn, int *nbrs, double *alpha,
	double *weights, int *scheme, int *inter, int *closest,	double *newline);

void mmult(double *x, double *y, int *nrx, int *ncx, int *ncy, double *ans);

void mycbind(double *a, double *b, int *ra, int *ca, int *cb, double *c);

void mycpyd(double *a, int *len, double *b);
void mycpyi(int *a, int *len, int *b);

void mydiag(double *d, int *n, int *ones, double *m);

void mymatchd(double *numa, double *numb, int *lnuma, int *lnumb, int *pos);
void mymatchi(int *numa, int *numb, int *lnuma, int *lnumb, int *pos);

void mymaxd(double *num, int *lnum, double *max, int *pos);
void mymaxi(int *num, int *lnum, int *max, int *pos);

void mymind(double *num, int *lnum, double *min, int *pos);
void mymini(int *num, int *lnum, int *min, int *pos);

void myrbind(double *a, double *b, int *ra, int *ca, int *rb, double *c);

void myrevd(double *dx, int *n, double *dy);
void myrevi(int *a, int *la, int *b);

void mysortd(double *a, int *la, double *sorted, int *order, int *inc);
void mysorti2(int *a, int *la, int *sorted, int *order, int *inc);

void mysvd(double *a, int *n, double *rvalues, double *rvectors,
  int *decreasing);

void myt(double *a, int *ra, int *ca, double *ta);

void mywhichd(double *num, int *lnum, double *a, int *pos);
void mywhichi(int *num, int *lnum, int *a, int *pos);

void pointsupdate(double *X, double *coeff, int *nn, int *index, int *remove,
	 int *pointsin, double *wts, double *l, int *N, double *alpha, int *r);

void pts(double *input, double *start, int *n, double *X);

void quadpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove,
	int *inter, int *nn, double *weights, int *docoeff);

void rmatsolve(double *m, int *n, double *inv);

/**
 * takes in the fwtnp lifting coefficient array, of size length(removelist)
 * by 3*max(n_r)+5 (zero filled,remove,nn,nbrs,alpha,weights,scheme,int,closest)
 *
 *    - matno is length(rem), nc is ncol(lca)
 **/
void transmatdual(double *lca, int *po, int *matno, int *lpo, int *nc,
	double *W, int *re);

void undopointsupdate(double *coeff, int *nbrs, int *index, int *remove, int *r,
	int *N, double *gamweights, double *l, double *lr, double *alpha, int *nn);

/**
 * adds newline to lca, adding zeros where appropriate.
 * NOTE: ncol(lca) is (3*nmax+5).
 * initial nr is nrow(lca)
 * initial nmax is ncol(lca)
 **/
void updatelca(double *lca, int *nr, int *nc, double *newline, double *newlca);

/**
 * gets schemehist from lca:
 *    1: LP
 *    2: QP
 *    3: CP
 **/
void schfromlca(double *lca, int *nr, int *nc, int *sch);

/**
 * gets interhist from lca:
 *    0: no intercept
 *    1: intercept
 **/
void interfromlca(double *lca, int *nr, int *nc, int *inter);

/**
 * gets clohist from lca:
 *    0: closest=FALSE
 *    1: closest=TRUE
 **/
void clofromlca(double *lca, int *nr, int *nc, int *clo);

/**
 * row no is given as R index (C index +1)
 **/
void nbrsfromlca(double *lca, int *nc, int *rowno, int *nbrs);

/**
 * row no is given as R index (C index +1)
 **/
void afromlca(double *lca, int *nc, int *rowno, double *alpha);

/**
 * row no is given as R index (C index +1)
 **/
void wfromlca(double *lca, int *nc, int *rowno, double *weights);

void invtnp(double *X, double *coeff, double *lengths, double *lengthsremove,
	int *pointsin, double *lca, int *nadd, int *N, int *lr, int *nc, int *outpo,
	double *outlen);

void findadds(int *rem, int *l, double *lca, int *nc, int *index, int *li,
	int *a);

void delrow(double *M, int *nr, int *nc, int *i, double *Mnew);

void getnbrs(double *X, int *remove, int *pointsin, int *lpo, int *neigh,
	int *closest, int *nbrs, int *index, int *nn);

/* this is a brute force, slightly inefficient way of
 sorting integers with index since nothing else seems to work */
void mysorti(int *a, int *la, int *sorted, int *order, int *inc);

void mmult2(double *A, double *B, int *ra, int *ca, int *cb, double *C);
void mmult3(double *A, double *B, int *ra, int *ca, int *cb, double *C);

#endif
