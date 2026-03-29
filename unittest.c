
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
#include<unistd.h>
#include<math.h>
#include<time.h>
#include<string.h>

#include "array.h"
#include "overlap.h"
#include "minunit.h"


int getopt(int argc, char *const argv[], const char *optstring);
extern char *optarg;
extern int optind, optopt;



int main(int argc, char *argv[])
{

	int k, sum;
	int p, K, resN, method, lim, fail, sph, g, nf;
	double emax, PiLow, Ubound;
	double BarOmega, MaxOmega;

	int *rcMax, *Nk;
	double *Pi, *pars;
	double **Mu, **OmegaMap, **Y;
	double ***S;
	
	double tol;
	int cond, cond1, cond2;
	
	MAKE_VECTOR(pars, 2);


/* DEFAULT PARAMETERS */

	p = 5;
	K = 7;
	
	method = 0; /* 0 - average; 1 - maximum */

	BarOmega = 0.0;
	MaxOmega = 0.0;
	sph = 0; /* 0 - nonspherical, 1 - spherical */
	emax = 0.90; /* Max eccentricity */
	PiLow = 1.0; /* Smallest mixing proportion */
	Ubound = 1.0; /* Means are distributed by U(0,Ubound) */
	resN = 100; /* Max # of dataset resimulations */
	pars[0] = 0.0000001; /* eps */
	pars[1] = 0.0000001; /* acc (Davies) */
	lim = 1000000; /* lim (Davies) */
	g = 100; /* sample size for generated datasets */
	nf = 1; /* # of simulated mixtures */

	char PIfname[100] = "TEST/Pi.test";
	char MUfname[100] = "TEST/Mu.test";
	char Sfname[100] = "TEST/LTSigma.test";
				
	
	srand(time(NULL));

	MAKE_VECTOR(Pi, K);
	MAKE_MATRIX(Mu, K, p);
	MAKE_3ARRAY(S, K, p, p);

	MAKE_MATRIX(OmegaMap, K, K);
	MAKE_VECTOR(rcMax, 2);
	
	MAKE_MATRIX(Y, g, p);
	MAKE_VECTOR(Nk, K);
	
/* RUNNING A CORRESPONDING FUNCTION */

	tol = 0.000001;


/* FUNCTION freadParameters */

	freadParameters(p, K, Pi, Mu, S, PIfname, MUfname, Sfname, 1);
	
	cond = (p == 5) && (K == 7) && (fabs(Pi[0] - 0.1831131) < tol) && (fabs(Pi[6] - 0.1569579) < tol) &&
			(fabs(Mu[0][0] - 0.0012522) < tol) && (fabs(Mu[1][1] - 0.6093170) < tol) &&
			(fabs(Mu[2][2] - 0.3875023) < tol) && (fabs(Mu[3][3] - 0.9270482) < tol) &&
			(fabs(Mu[4][4] - 0.7657742) < tol) && (fabs(Mu[6][4] - 0.7767524) < tol) &&
			(fabs(S[0][0][0] - 0.0162117) < tol) && (fabs(S[1][0][1] + 0.0005859) < tol) &&
			(fabs(S[2][0][2] + 0.0060184) < tol) && (fabs(S[4][3][0] - 0.0047450) < tol) &&
			(fabs(S[5][0][4] - 0.0008476) < tol) && (fabs(S[6][4][4] - 0.0208782) < tol);
	printf("function \"freadParameters\"   -  ");
	func_test(cond, "OK\n", "Test failed\n");


/* FUNCTION ExactOverlap */
		
	ExactOverlap(p, K, Pi, Mu, S, pars, lim, OmegaMap, &BarOmega, &MaxOmega, rcMax);

	cond = (fabs(BarOmega - 0.017414) < tol) && (fabs(MaxOmega - 0.198721) < tol);
	printf("function \"ExactOverlap\"      -  ");
	func_test(cond, "OK\n", "Test failed\n");


/* FUNCTION genData */

	Anull(Y, g, p);
	genData(p, K, Pi, Mu, S, g, Y, Nk);
	
	sum = 0;
	for (k=0; k < K; k++){
		sum = sum + Nk[k];
	}
	
	cond = (g == sum) && (Y[0][0] != 0.0) && (Y[g-1][p-1] != 0.0);
	
	printf("function \"genData\"           -  ");
	func_test(cond, "OK\n", "Test failed\n");


/* FUNCTION OmegaClust */
	
	// test 1
	sph = 0;
	method = 1;
	OmegaClust(0.15, method, p, K, PiLow, Ubound, emax, pars, lim, resN, sph, Pi, Mu, S, OmegaMap, &BarOmega, &MaxOmega, rcMax, &fail);
	cond1 = (fabs(0.15 - MaxOmega) < tol);

	// test 2
	PiLow = 0.1;
	sph = 1;
	method = 0;
	OmegaClust(0.03, method, p, K, PiLow, Ubound, emax, pars, lim, resN, sph, Pi, Mu, S, OmegaMap, &BarOmega, &MaxOmega, rcMax, &fail);
	cond2 = (fabs(0.03 - BarOmega) < tol) && (fabs(S[0][1][0]) < tol) && (fabs(S[1][2][1]) < tol);
	
	cond = (cond1 && cond2);
	printf("function \"OmegaClust\"        -  ");
	func_test(cond, "OK\n", "Test failed\n");


/* FUNCTION OmegaClust */

	sph = 0;
	BarOmega = 0.03;
	MaxOmega = 0.15;
	OmegaBarOmegaMax(p, K, PiLow, Ubound, emax, pars, lim, resN, sph, Pi, Mu, S, OmegaMap, &BarOmega, &MaxOmega, rcMax, &fail);

	cond = (fabs(BarOmega - 0.03) < tol) && (fabs(MaxOmega - 0.15) < tol);
	printf("function \"OmegaBarOmegaMax\"  -  ");
	func_test(cond, "OK\n", "Test failed\n");


	FREE_MATRIX(Y);
	FREE_VECTOR(Nk);

	FREE_VECTOR(Pi);
	FREE_MATRIX(Mu);
	FREE_3ARRAY(S);

	FREE_MATRIX(OmegaMap);
	FREE_VECTOR(rcMax);

	FREE_VECTOR(pars);

	return 0;

}

