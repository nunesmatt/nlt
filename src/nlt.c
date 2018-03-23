#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <Rinternals.h>

/*
  TO-DO:
    Rewrite:
      intervals
      mysortd

  Optimisation hints:
    rather than for(){dofunct()} write dofunc{for{}}
    Pass structures by reference, not by value. and reduce the number of function parameters
    Use shift operations >> and << instead of integer multiplication and division, where possible.
    If you repeatedly divide by x, consider computing 1/x and multiplying by the result

*/

#include "adlift.h"
#include "nlt.h"

#define MULTIPLY2(a) a << 1
#define MULTIPLY3(a) a + (a << 1)
#define MULTIPLY4(a) a << 2
#define DIVIDE2(a) a >> 1
#define DIVIDE4(a) a >> 2

struct theConstants {
  double *input;
  double *f;
  int nkeep;
  int intercept;
  int initboundhandl;
  int neighbours;
  int closest;
  int LocalPred;
  int n;
  int *traj;
  bool doW;
  bool varonly;
  bool do_w_or_v;
} consts;

struct theOutputs {
  double *coeff; //output
  double *lengthsremove; //output
  double *lengths; //output
  int *pointsin; //output
  double *W; //output
  double *v; //output
  double *lca; //DONT KNOW WHAT THIS IS FOR
  int *nc; //DONT KNOW WHAT THIS IS FOR
} outputs;

void calculateIntervals(double *X){
  double I[consts.n + 1];
  intervals(X,&consts.initboundhandl,&consts.n,I);
  for(int i=0; i<consts.n; i++){
    outputs.lengths[i] = I[i+1] - I[i];
  }
}

void calculateCoeff(double *X){
  double sX[consts.n];
  int one = 1;
  mysortd(X, &consts.n, sX, outputs.pointsin, &one);
  memcpy(outputs.coeff, consts.f, consts.n * sizeof(double));
}

struct data{
  int scheme;
  int nn;
  int r;
  int nr;
  int *nbrs;
  int *index;
  double *weights;
  double *alpha;
};

struct data findNeighbours(double *X, int *remove, int *N){
  //int nnmax = 2 * consts.neighbours;
  int nnmax = MULTIPLY2(consts.neighbours);
  struct data d1;
  d1.nbrs = (int *)malloc(nnmax * sizeof(int));
  d1.index = (int *)malloc(nnmax * sizeof(int));
  d1.nn = 0;
  getnbrs(X, remove, outputs.pointsin, N, &consts.neighbours,
    &consts.closest, d1.nbrs, d1.index, &d1.nn);

  struct data d2;
  d2.nbrs = (int *)malloc(d1.nn * sizeof(int));
  d2.index = (int *)malloc(d1.nn * sizeof(int));
  d2.weights = (double *)malloc(d1.nn * sizeof(double));
  d2.alpha = (double *)malloc(d1.nn * sizeof(double));
  d2.nn = d1.nn;
  memcpy(d2.nbrs, d1.nbrs, d2.nn * sizeof(int));
  memcpy(d2.index, d1.index, d2.nn * sizeof(int));

  free(d1.nbrs);
  free(d1.index);

  return d2;
}

struct data fitCurve_1_to_4(double *X, int *remove, int *N){
  struct data d = findNeighbours(X, remove, N);
  int one = 1;
  d.scheme = consts.LocalPred;
  switch(consts.LocalPred){
    case 1:
      linearpred(outputs.pointsin, X, outputs.coeff, d.nbrs, remove,
        &consts.intercept, &d.nn, d.weights, &one);
      break;
    case 2:
      quadpred(outputs.pointsin, X, outputs.coeff, d.nbrs, remove,
        &consts.intercept, &d.nn, d.weights, &one);
      break;
    case 3:
      cubicpred(outputs.pointsin, X, outputs.coeff, d.nbrs, remove,
        &consts.intercept, &d.nn, d.weights, &one);
      break;
    case 4:
      d.scheme=1;
      adaptpred(outputs.pointsin, X, outputs.coeff, d.nbrs, remove,
        &consts.intercept, &d.nn, d.weights, &d.scheme, &one);
      break;
  }
  return d;
}

struct data fitCurve_5(double *X, int *remove, int *N){
  //int nnmax = 2 * consts.neighbours;
  int nnmax = MULTIPLY2(consts.neighbours);
  int one = 1;
  struct data d;
  d.nbrs = (int *)malloc(nnmax * sizeof(int));
  d.index = (int *)malloc(nnmax * sizeof(int));
  d.weights = (double *)malloc(nnmax * sizeof(double));
  d.nn = 0;
  d.scheme=1;
  adaptneigh(outputs.pointsin, X, outputs.coeff, d.nbrs, remove,
    &consts.intercept, &d.nn, d.weights, &d.scheme, &consts.closest, d.index,
    &consts.neighbours, N, &one);
  return d;
}

double * lca(struct data d, int j, int *remove){
  //double newline[3 * d.nn + 5];
  double newline[MULTIPLY3(d.nn) + 5];
  makelcaline(remove, &d.nn, d.nbrs, d.alpha, d.weights, &d.scheme,
    &consts.intercept, &consts.closest, newline);
  //int max = (*outputs.nc>=(3*d.nn+5)) ? *outputs.nc : (3*d.nn+5);
  int max = (*outputs.nc>=(MULTIPLY3(d.nn)+5)) ? *outputs.nc : (MULTIPLY3(d.nn)+5);
  double *tmplca = (double *)malloc(max * j * sizeof(double));
  d.nr = 0;
  updatelca(outputs.lca, &d.nr, outputs.nc, newline, tmplca);
  return tmplca;
}

double * makeWnew(double *Wnew, int j, int r){
  int dim = consts.n - j + 1;
  int dim1 = dim - 1;
  int dim2 = dim * consts.n;
  double *Wtmp = (double *)malloc(dim2 * sizeof(double));
  memcpy(Wtmp, Wnew, dim2);
  free(Wnew);
  double *Wrtn = (double *)malloc(dim1 * consts.n * sizeof(double));
  delrow(Wtmp,&dim,&consts.n,&r,Wrtn);
  free(Wtmp);
  return Wrtn;
}

void updateLenghtsAndPoints(int *N, int r){
  double len2[*N];
  int po[*N];
  memcpy(len2, outputs.lengths, *N * sizeof(double));
  memcpy(po, outputs.pointsin, *N * sizeof(int));
  getridd(len2,N,&r,outputs.lengths);
  getridi(po,N,&r,outputs.pointsin);
}

/**
 *
 * This function is called 2010000 times per point in each table in Kathryns code. 11 * 8 points per simulation = 176,880,000.
 *
 * _input: location of sample
 * _f: gain at specified location
 * _nkeep: The number of scaling coefficients to be kept in the final
           representation of the initial signal. This must be at least
           two.
 * _intercept: Indicates whether or not the regression curve includes an
               intercept.
 * _initboundhandl: variable specifying how to handle the boundary at the
          start of the transform.  Possible values are ‘"reflect"’ -
          the intervals corresponding to the first and last datapoints
          are taken to have the respective grid values as midpoints;
          and ‘"stop"’ - the first and last intervals have the first
          and last grid values (respectively) as outer endpoints.
 * _neighbours: The number of neighbours over which the regression is
          performed at each step. If closest is false, then this in
          fact denotes the number of neighbours on each side of the
          removed point.
 * _closest: Refers to the configuration of the chosen neighbours. If
          ‘closest’ is false, the neighbours will be chosen
          symmetrically around the removed point. Otherwise, the
          closest neighbours will be chosen.
 * _LocalPred: The type of regression to be performed. Possible options are
          ‘LinearPred’, ‘QuadPred’, ‘CubicPred’, ‘AdaptPred’ and
          ‘AdaptNeigh’
 * _n: lengthof(_input)
 * _coeff: Main output1 (of length _n): vectorvector of detail and scaling coefficients in the wavelet
           decomposition of the signal.
 * _lengthsremove: Main output2 (of length _n - _nkeep): vector of interval lengths corresponding to the points
          removed during the transform (in ‘removelist’).
 * _lengths: vector of (updated) interval lengths at the end of the
          transform. This is of length ‘nkeep’
 * WHAT IS THIS? _lca: 2D array with _n-_nkeep rows, and 6 * _neighbours + 5 cols
 * _pointsin: OUTPUT Vector of length _nkeep. Indices into ‘X’ of the scaling coefficients in the wavelet
          decomposition. These are the indices of the ‘X’ values which
          remain after all points in ‘removelist’ have been predicted
          and removed. This has length ‘nkeep’.
 * _nc: WHAT IS THIS? is always 0. when called from R is called with as.integer(0)
 * _traj: Vector of length (length(‘x’)-‘keep’). It gives the
          trajectory for the modified lifting algorithm to follow, i.e.
          it gives the order of point removal. (mod in R)
 * _doW: A boolean indicating whether the transform matrix should be
          computed and returned. Combined with varonly to create ex
 * _W: 2D array of doubles. Dimensions _n * _n. The transform matrix
 * _varonly: A boolean indicating whether only the coefficient variances
          should be returned (if ‘do.W=TRUE’). Combines with _doW to create ex
 * _v: vector of type double of length _n. The coefficient variances
 **/
void fwtnpperm(double *_input,
               double *_f,
               int *_nkeep,
               int *_intercept,
               int *_initboundhandl,
               int *_neighbours,
               int *_closest,
               int *_LocalPred,
               int *_n,
               double *_coeff,
               double *_lengthsremove,
               double *_lengths,
               double *_lca,
               int *_pointsin,
               int *_nc,
               int *_traj,
               int *_doW,
               double *_W,
               int *_varonly,
               double *_v)
{
  consts.input = _input;
  consts.f = _f;
  consts.nkeep = *_nkeep;
  consts.intercept = *_intercept;
  consts.initboundhandl = *_initboundhandl;
  consts.neighbours = *_neighbours;
  consts.closest = *_closest;
  consts.LocalPred = *_LocalPred;
  consts.n = *_n;
  consts.traj = _traj;
  consts.doW = (bool)*_doW;
  consts.varonly = (bool)*_varonly;
  consts.do_w_or_v = (bool)(*_doW + *_varonly);

  if(*_doW + *_varonly > 1){
    printf("Cannot have both doW and varonly set.\n Exiting");
    exit(1);
  }

  outputs.coeff = _coeff;
  outputs.lengthsremove = _lengthsremove;
  outputs.lengths = _lengths;
  outputs.pointsin = _pointsin;
  outputs.W = _W;
  outputs.v = _v;
  outputs.lca = _lca;
  outputs.nc = _nc;

  double X[consts.n];
  memcpy(X, consts.input, consts.n * sizeof(double));
  calculateIntervals(X);
  calculateCoeff(X);

  double *Wnew = (double *)malloc(consts.n * consts.n * sizeof(double));
  if(consts.do_w_or_v){
    for(int i=0 ; i < consts.n; i++){
      Wnew[(i*consts.n)+i]=1;
    }
  }

  int N = consts.n;
  if(consts.nkeep != consts.n){
    for (int j=1; j <= (consts.n-consts.nkeep); j++) {
      int remove = consts.traj[j-1];
      struct data d;
      if(consts.LocalPred == 5){
        d = fitCurve_5(X, &remove, &N); //nn doesnt get initialised!
      } else {
        d = fitCurve_1_to_4(X, &remove, &N);
      }
      pointsupdate(X,outputs.coeff, &d.nn, d.index, &remove, outputs.pointsin,
        d.weights, outputs.lengths, &N, d.alpha, &d.r);
      outputs.lengthsremove[j-1]= outputs.lengths[d.r-1];

      double * tmplca = lca(d, j, &remove);

      if(consts.do_w_or_v){
        if(consts.varonly){
          for(int i=0; i < consts.n; i++){
            for(int k=0; k < d.nn; k++){
              Wnew[(d.r-1)*consts.n+i]-=
                d.weights[k] * Wnew[(d.index[k]-1)*consts.n+i];
            }
          }
          for(int i=0; i < consts.n; i++){
            for(int k=0; k < d.nn; k++){
              Wnew[(d.index[k]-1)*consts.n+i] +=
                d.alpha[k] * Wnew[(d.r-1)*consts.n+i];
            }
            outputs.v[remove-1] += pow(Wnew[(d.r-1)*consts.n+i],2);
          }
          Wnew = makeWnew(Wnew, j, d.r);
        }
        else {
          for(int i = 0; i < consts.n; i++){
            for(int k = 0; k < d.nn; k++){
              Wnew[(remove-1) * consts.n + i] -=
                d.weights[k] * Wnew[(d.nbrs[k]-1) * consts.n + i];
            }
          }
          for(int i=0; i < consts.n; i++){
            for(int k=0; k< d.nn; k++){
              Wnew[(d.nbrs[k]-1) * consts.n + i] +=
                d.alpha[k] * Wnew[(remove-1) * consts.n + i];
            }
          }
        }
      }


      free(d.nbrs);
      free(d.alpha);
      free(d.weights);
      free(d.index);

      updateLenghtsAndPoints(&N, d.r);

      int dim= d.nr * *outputs.nc;
      memcpy(outputs.lca, tmplca, dim * sizeof(double));
      free(tmplca);
      N-=1;
    } //end of j loop
  } //end of if(consts.nkeep != consts.n){

  if(consts.do_w_or_v){
    if(consts.varonly){
      for(int i = 0; i < N; i++){
        for(int k = 0; k < consts.n; k++){
          outputs.v[outputs.pointsin[i]-1] += pow(Wnew[(i*consts.n)+k],2);
        }
      }
      int dim1= consts.nkeep * consts.n;
      memcpy(outputs.W, Wnew, dim1 * sizeof(double));
    }
    else{
      int dim2sq = pow(consts.n,2);
      memcpy(outputs.W, Wnew, dim2sq * sizeof(double));
    }
    free(Wnew);
  }

  //these got assigned new memory????
  _W = outputs.W;
  _v = outputs.v;
}

//memcpy(dst, src, no.Bytes)

/*int main(int argc, char const *argv[]) {
  printf("MAIN");
  return 0;
}*/
