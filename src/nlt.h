#ifndef __ADLIFT
#define __ADLIFT

void fwtnpperm(double *input,
               double *f,
               int *nkeep,
               int *intercept,
               int *initboundhandl,
               int *neighbours,
               int *closest,
               int *LocalPred,
               int *n,
               double *coeff,
               double *lengthsremove,
               double *lengths,
               double *lca,
               int *pointsin,
               int *nc,
               int *traj,
               int *doW,
               double *W,
               int *varonly,
               double *v);

#endif
