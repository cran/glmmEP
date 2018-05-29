/*
 *  Part of R package glmmEP
 *  Copyright (C) 2018  M.P. Wand and J.C.F. Yu
 *
 *  Unlimited use and distribution.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void F77_SUB(asn)(double *a1, double *a2, double *A2ina1, int *idmn,
                  int *idmnsq, int *lena2, double *A2mat,
                  double *xm2A2, double *DdPlus, double *wkv,
                  int *ipvt, double *det, double *work, double *ans);

void F77_SUB(cpbt)(double *a1, double *a2, double *b1, double *b2,
                   double *c0, double *c1, int *idmn, int *idmsq,
                   int *lena2, double *Dd, double * DdPlus, 
                   double *wka, double *wkb, double *A2ina1, 
                   double *B2inb1, double *A2inc1, int *ipvt, 
                   double *det, double *work, double *A2mat,
                   double *A2neg, double *B2mat, double *B2neg, 
                   double *ans);

void F77_SUB(epllk)(double *beta,double *etaSg0,double *etaSg,
                    int *m,int *nVec, int *nMax, int *numObs,
                    int *indStt,int *idF,int *idR, int *idRsq,
                    int *lena, int *lena2, int*nlena, double *yDagg,
                    double *Xf, double *Xr, double *yDCur, 
                    double *XfCur, double *XrCur, double *uHat,
                    double *Dd, double *DdPlus, double *Xbeta,
                    double *XuHat, double *etaFtS, double *etaStF,
                    double *SUMlt, double *c1Cur, double *etaIN1,
                    double *etaIN2, double *eINa1, double *eINa2,
                    double *eINb1, double *eINb2, double *etaPvF,
                    double *etaPvS,double *etaCrF, double *etaCrS, 
                    double *wk1, int *ipvt, double *A2ina1,
                    double *A2inc1, double *A2mat, double *A2str,
                    double *R2comp, double *wk2,double *R5, double *R5TA2, 
                    double *vR5TA2, double *xkpan1, double *xkpan2, 
                    double *B2inb1, double *work, double *A2neg, 
                    double *B2mat, double *B2str, double *B2neg, 
                    double *xm2A2, double *det, double *wka, double *wkb, 
                    double *wkv, double *xmiscl, double *etaOut);

void F77_SUB(kpbt)(double *a1, double *a2, double *c0, double *c1,
                   int *idmn, int *idmsq, int *lena2, double *Dd,
                   double *DdPlus, double *wk1, double *A2ina1,
                   double *A2inc1, int *ipvt, double *A2mat, 
                   double *A2str, double *R2comp, double *wk2,
                   double *R5, double *R5TA2, double *vR5TA2,
                   double *ans1, double *ans2);


void F77_SUB(logphi)(double *x, double *ans);

void F77_SUB(logdet)(double *A, int *idmn, int *ipvt, double *work,
                     double *det, double *ans);

void F77_SUB(zetad)(double *x, double *ans);

static const R_FortranMethodDef FortEntries[] = {
    {"asn",    (DL_FUNC) &F77_SUB(asn),   14},
    {"cpbt",   (DL_FUNC) &F77_SUB(cpbt),  24},
    {"epllk",  (DL_FUNC) &F77_SUB(epllk), 65},
    {"kpbt",   (DL_FUNC) &F77_SUB(kpbt),  22},
    {"logphi", (DL_FUNC) &F77_SUB(logphi), 2},
    {"logdet", (DL_FUNC) &F77_SUB(logdet), 6},
    {"zetad",  (DL_FUNC) &F77_SUB(zetad),  2},
    {NULL, NULL, 0}
};

void R_init_glmmEP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
