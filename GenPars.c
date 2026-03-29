
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





#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"

#include "overlap.h"


/* genSigma : generates covariance matrix based on (p + 1) observations (unstable covariance matrix)
 * Parameters:
 * 		p - number of dimensions
 * 		VC - variance-covariance matrix
 */

void genSigma(int p, double **VC){

	int i,j,k,n;
	double **x, *mu;

	n = p + 1;

	MAKE_MATRIX(x, n, p);
	MAKE_VECTOR(mu, p);

	anull(mu, p);

	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			x[i][j] = rnor(0.0, 1.0);
			mu[j] = mu[j] + x[i][j];
		}
        }

	for (j=0;j<p;j++){
		mu[j] = mu[j] / n;
	}	
	
	Anull(VC, p, p);

	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			for (k=0; k<p; k++){
				VC[j][k] = VC[j][k] + (x[i][j] - mu[j]) * (x[i][k] - mu[k]);
			}
		}
        }

	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			VC[j][k] = VC[j][k] / (n - 1);
		}
	}

	FREE_MATRIX(x);
	FREE_VECTOR(mu);

}


/* genSigmaEcc : generates covariance matrix  with prespecified eccentricity
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		emax - maximum eccentricity
 * 		S - set of variance-covariance matrices
 */

void genSigmaEcc(int p, int K, double emax, double ***S){

	int i, k;

	double dtmt, minL, maxL, e;
	double *Eig;
	double **VC, **L, **R;

	MAKE_VECTOR(Eig, p);
	MAKE_MATRIX(VC, p, p);
	MAKE_MATRIX(L, p, p);
	MAKE_MATRIX(R, p, p);

	for (k=0; k<K; k++){
		genSigma(p, VC);
		cpy2(VC, p, p, S, k);		

		cephes_symmeigens_down(p, Eig, VC, &dtmt);

/*		for (i = 0; i < p; i++) printf("%lf ", Eig[i]);
		printf("\n"); */


		i = vecMin(Eig, p, &minL);
		
/*		printf("mini = %d\n", i); */

		i = vecMax(Eig, p, &maxL);

/*		printf("maxi = %d\n", i);*/

		e = pow(1 - minL / maxL, 0.5);


		if (e > emax){

			Anull(L, p, p);

			for (i=0; i<p; i++){
				Eig[i] = maxL * (1 - emax * emax * (maxL - Eig[i]) / (maxL - minL));
				L[i][i] = Eig[i];
			}

			XAXt(VC, p, L, R);
			cpy2(R, p, p, S, k);

		}

	}

	FREE_MATRIX(VC);
	FREE_MATRIX(L);
	FREE_MATRIX(R);

}



/* genSphSigma : generates spherical covariance matrix
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		S - set of variance-covariance matrices
 */

void genSphSigma(int p, int K, double ***S){

	int i, k;
	double r;
	double **L;

	MAKE_MATRIX(L, p, p);		

	Anull(L, p, p);
	
	for (k=0; k<K; k++){

		r = runir(0.0, 1.0);
		for (i=0; i<p; i++){			
			L[i][i] = r;
		}

		cpy2(L, p, p, S, k);
	
	}

	FREE_MATRIX(L);

}



/* genSphSigma : generates matrix of means
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		Mu - set of mean vectors
 * 		Ubound - upper bound for the hypercube
 */

void genMu(int p, int K, double **Mu, double Ubound){
		
	int i, k;
	
	if (Ubound <= 0) Ubound = 1.0;
	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
		
			Mu[k][i] = runir(0.0, Ubound);

		}
	}

}



/* genPi : generates mixing proportions
 * Parameters:
 * 		K - number of components
 * 		PiLow - smallest possible mixing proportion
 * 		Pi - vector of mixing proportions
 */

void genPi(int K, double PiLow, double *Pi){

	int flag, k;
	double s;

	flag = 0;


	if ((PiLow >= 1) | (PiLow <= 0)){
/*		printf("Warning: PiLow is out of range... generated equal mixing proportions...\n"); */
		for (k=0; k<K; k++){
			Pi[k] = 1.0 / K;
		}
	} else {
		s = 0.0;
		for (k=0; k<K; k++){
			Pi[k] = rgamma(1.0);
			s += Pi[k];
		}
		for (k=0; k<K; k++){
			Pi[k] = PiLow + Pi[k] / s * (1 - K * PiLow);
			if (Pi[k] < PiLow){
				flag = 1;
				break;
			}
		}
		if (flag == 1){
/*			printf("Warning: PiLow is too high... generated equal mixing proportions...\n"); */
			for (k=0; k<K; k++){
				Pi[k] = 1.0 / K;
			}
		}
		
		
	}

}
