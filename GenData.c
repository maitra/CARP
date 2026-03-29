
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


/* genData: generates a dataset given mixture parameters
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		Pi - vector of mixing proportions
 * 		Mu - set of mean vectors
 * 		S - set of covariance matrices
 * 		n - number of observations
 * 		Y - simulated dataset
 * 		Nk - cluster sizes
 */

void genData(int p, int K, double *Pi, double **Mu, double ***S, int n, double **Y, int *Nk){

	int i, j, k;
	int curr;
	double **X;

	MAKE_MATRIX(X, n, p);

	rmulti(n, K, Pi, Nk);

	curr = 0;
	for (k=0; k<K; k++){	
		
		rMVN(Nk[k], p, Mu[k], S[k], X);
		
		for (i=0; i<Nk[k]; i++){
			for (j=0; j<p; j++){
				Y[curr+i][j] = X[i][j];
			}
		}

		curr = curr + Nk[k];

	}

	FREE_MATRIX(X);

}

