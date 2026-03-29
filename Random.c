
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




#include <stdio.h>
#include <stdlib.h>	
#include <math.h>
#include "array.h"
#include "overlap.h"

#define PI 3.141592653589793
#define EE 2.718281828459046

long int random(void);   /* declare in order to avoid undef with ISO */
void srandom(unsigned int seed); /* declare in order to avoid undef with ISO */
double rxxx1;
int ixxx1=0;
/*************************************/
void setseed(unsigned int s)
{
  srandom(s);
}
/*************************************/
long genseed(void)
{
  return((long)random());
}
/*************************************/
double runi(void)
{
  return ((double)random() + 0.5)/((double)RAND_MAX + 1.0);
}
/*************************************/
double runir(double a,double b)
{
  return (b-a)*runi()+a;
}
/*************************************/
double rnor(double mu,double sd)
{
   double e,v1,v2,w;
   if(ixxx1==0){
      do{
         v1=2*runi()-1.;
         v2=2*runi()-1.;
         w=v1*v1+v2*v2;      
      }while(w>1.);
      e=sqrt((-2.*log(w))/w);
      rxxx1=v1*e;
      ixxx1=1;
      return v2*e*sd+mu;
   }
   else{
      ixxx1=0;
      return rxxx1*sd+mu;
   }
}
/*************************************/
double rexp(double lambda)
{
  return -log(runi())/lambda;
}
/*************************************/
double rgamma(double alpha)
{
  double r1,r2,aa,x,w,c1,c2,c3,c4,c5;
  if(alpha<=0.)return 0.;
  if(alpha == 1.)return rexp(1.);
  if(alpha<1){
    aa=(alpha+EE)/EE;
    do{
      r1=runi();
      r2=runi();
      if(r1>1./aa){
	x = -log(aa*(1.-r1)/alpha);
	if(r2<pow(x,(alpha-1.)))return x;
      }
      else{
	x = pow((aa*r1),(1./alpha));
	if(r2<exp(-x))return x;
      }
    }while(r2<2);
   }
  else{
    c1=alpha-1;
    c2=(alpha-1./(6.*alpha))/c1;
    c3=2./c1;
    c4=c3+2.;
    c5=1./sqrt(alpha);
    do{
      do{
	r1=runi();
	r2=runi();
	if(alpha>2.5)r1=r2+c5*(1.-1.86*r1);
      }while(r1<=0 || r1 >= 1);
      w=c2*r2/r1;
      if(c3*r1+w+1/w <= c4)return c1*w;
      if(c3*log(r1)-log(w)+w<1)return c1*w;
    }while(r2<2);
  }
  exit(1);
}
/*************************************/
/* Generation of random realizations from Multinomial distribution
 *	n  - number of trials 
 *  K  - dimensionality
 *  pi - multinomial proportions
 *  Nk - produced random realization
 */
void rmulti(int n, int K, double *pi, int *Nk){
   
   int i, k;
   double X;
   double *cumPI;
     

   MAKE_VECTOR(cumPI, K);

   for(k=0; k < K; k++){
	   Nk[k] = 0;
   }

   cumPI[0] = pi[0];
   for(k=1; k < K; k++){
	   cumPI[k] = cumPI[k-1] + pi[k];
   }

   for (i = 0; i < n; i++){
	   X = runi();
	   for(k=0; k < K; k++){
		   if (X < cumPI[k]){
			   Nk[k] = Nk[k] + 1;
			   break;
		   }
	   }
   }

   FREE_VECTOR(cumPI);

}
/*************************************/
/* Generation of random realizations from multivariate Gaussian distribution
 *	n  - sample size
 *  K  - dimensionality
 *  MU - mean vector
 *  S  - covariance matrix
 *  X  - produced random realizations
 */
void rMVN(int n, int p, double *MU, double **S, double **X){

	int i, j;
	double dtmt;

	double *Eig;
	double **Ga, **L, **Sh, **Z;

	MAKE_MATRIX(Sh, p, p);
	MAKE_VECTOR(Eig, p);
	MAKE_MATRIX(Ga, p, p);	
	MAKE_MATRIX(L, p, p);
	MAKE_MATRIX(Z, n, p);


	cpy(S, p, p, Ga);
	cephes_symmeigens_down(p, Eig, Ga, &dtmt);

	Anull(L, p, p);
	for (i=0; i<p; i++){
		L[i][i] = pow(Eig[i], 0.5);
	}
	XAXt(Ga, p, L, Sh);


	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			Z[i][j] = rnor(0.0, 1.0);
		}
	}

	multiply(Z, n, p, Sh, p, p, X);

	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			X[i][j] = X[i][j] + MU[j];
		}
	}

	FREE_MATRIX(Sh);
	FREE_MATRIX(Ga);
	FREE_MATRIX(L);
	FREE_VECTOR(Eig);
	FREE_MATRIX(Z);

}
