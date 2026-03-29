

/*

    CARP: Clustering Algorithms' Referee Package
    Copyright (C) 2010  Volodymyr Melnykov and Ranjan Maitra

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	Authors' contact information:
	Volodymyr Melnykov						Ranjan Maitra
	Volodymyr.Melnykov@ndsu.edu				maitra@iastate.edu
	Department of Statistics				Department of Statistics
	North Dakota State University			Iowa State University
	Fargo, ND 58102							Ames, IA 50011

*/




#ifndef OVERLAP_H
#define OVERLAP_H

void ExactOverlap(int p, int K, double *Pi, double **Mu, double ***S, double *pars, int lim, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax);
void OmegaClust(double Omega, int method, int p, int K, double PiLow, double Ubound, double emax, double *pars, int lim, int resN, int sph, double *Pi, double **Mu, double ***S, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax, int (*fail));
void OmegaBarOmegaMax(int p, int K, double PiLow, double Ubound, double emax, double *pars, int lim, int resN, int sph, double *Pi, double **Mu, double ***S, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax, int (*fail));
void dsyev_(char *JOBZp, char *UPLOp,int *Np, double *A, int *LDAp, double *Wp, double *WORK, int *LWORK, int *INFOp);
int vecMin(double *x, int p, double (*min));
int vecMax(double *x, int p, double (*max));
void Anull(double **X, int ax, int bx);
void anulli(int *x, int p);
void anull(double *X, int p);
void XAXt(double **X, int p, double **A, double **Res);
void XAXt2(double **X, int p, double **A, double ***Res, int k);
void cpy(double **a, int nrows, int ncols, double **b);
void cpy1(double ***a, int k, int nrows, int ncols, double **b);
void cpy2(double **a, int nrows, int ncols, double ***b, int k);
void tA(double **A, int a, int b, double **Res);
void matxvec(double **a, int arows, int acols, double *x, int xrows, double *y);
void multiply(double **a, int arows, int acols, double **b, int brows, int bcols, double **c);
void cephes_symmeigens_down(int p, double *eval, double **A, double (*determinant));
void cxS(int p, int K, double ***S, double c);
double qfc(double* lb1, double* nc1, int* n1, int *r1in, double *sigmain, double *c1in, int *lim1in, double *accin, double* trace, int* ifault);
void genPi(int K, double PiLow, double *Pi);
void genMu(int p, int K, double **Mu, double Ubound);
void genSigmaEcc(int p, int K, double emax, double ***S);
void genSphSigma(int p, int K, double ***S);
void fprintOverlap(int K, double **OmegaMap, double BarOmega, double MaxOmega, int *rcMax, int fileN, char *overmap, char *overbarmax);
void printOverlap(int K, double **OmegaMap, double BarOmega, double MaxOmega, int *rcMax);
void fprintParameters(int p, int K, double *Pi, double **Mu, double ***S, char *PIfname, char *MUfname, char *Sfname, int fileN);
void printParameters(int p, int K, double *Pi, double **Mu, double ***S);
void freadParameters(int p, int K, double *Pi, double **Mu, double ***S, char *PIfname, char *MUfname, char *Sfname, int fileN);
void fprintData(int n, int p, int K, double **x, int *Nk, int fileN, char *dataX, char *Nksizes);
void array1to2(int a, int b, double *y, double **x);
void array1to3(int a, int b, int c, double *y, double ***x);
void array2to1(int a, int b, double *y, double **x);
void array3to1(int a, int b, int c, double *y, double ***x);
void array3to1(int a, int b, int c, double *y, double ***x);
void setseed(unsigned int s);
long genseed(void);
double runir(double a,double b);
double rnor(double mu,double sd);
double rgamma(double alpha);
void rmulti(int n, int K, double *pi, int *Nk);
void rMVN(int n, int p, double *MU, double **S, double **X);
void genData(int p, int K, double *Pi, double **Mu, double ***S, int n, double **Y, int *id);

#endif /* OVERLAP_H */
