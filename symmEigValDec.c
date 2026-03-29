
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
#include <string.h>
#include "array.h"
#include "cephes_eigens.h"


void cephes_symmeigens_down(int p, double *eval, double **A, 
			    double (*determinant))

/*
  input:

  p - dimension of matrix
  A - pxp symmetric matrix (only lower triangle used) destroyed on return

  returns:
  eval - vector of eigenvalues (in ascending order)
  evec - matrix of eigenvectors 
  determinant - determinant of the symmetric matrix (as calculated by the 
                product of the eigenvalues

*/
	
{
	int i, j;
	double *As, *Evec, *Evalues;

	MAKE_VECTOR(As, p * (p + 1) / 2);

	for (i = 0; i < p; i++) {
		for (j = 0; j <= i; j++) As[(i * i + i)/2 + j] = A[i][j];
	}

	MAKE_VECTOR(Evec, p * p);	
	MAKE_VECTOR(Evalues, p);

	cephes_eigens(As, Evec, Evalues, p);
	
	for (i = 0; i < p; i++){
		eval[i] = Evalues[p - i - 1];
	}

	for (i = 0; i < p; i++) {
		for (j = 0; j < p; j++){
			A[j][p-i-1] = Evec[p * i + j];
		}
	}
	                        

/*
 	printf("EigenValues:\n");
	for (i = 0; i < p; i++){
		printf("EigVal: %lf \n", eval[i]);
	}
 	printf("EigenVectors:\n");
	for (i = 0; i < p; i++) {
		for (j = 0; j < p; j++){
			printf("%lf ", A[i][j]);
		}
		printf("\n");
	}
*/
	
	(*determinant)=1.0;
	
	for (i = 0; i < p; i++) (*determinant)*=eval[i];
	
	FREE_VECTOR(As);
	FREE_VECTOR(Evalues);
	FREE_VECTOR(Evec);

	return;
}

