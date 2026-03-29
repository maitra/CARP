
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

int getopt(int argc, char *const argv[], const char *optstring);
extern char *optarg;
extern int optind, optopt;




/*
	C-MixSim generates mixtures with prespecified level of overlap and datasets from these mixtures.
*/

int main(int argc, char *argv[])
{

	int i, k;
	int p, K, resN, method, lim, fail, sph, g, nf;
	double emax, PiLow, Ubound;
	double BarOmega, MaxOmega;

	int *rcMax;
	double *Pi, *pars;
	double **Mu, **OmegaMap;
	double ***S;
	
	int c, errflag;
	char **pp;
	
	pp = NULL;

	MAKE_VECTOR(pars, 2);

/* DEFAULT PARAMETERS */

	p = 2;
	K = 2;
	
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
	g = 0; /* sample size for generated datasets */
	nf = 1; /* # of simulated mixtures */

	char buff[200] = "";
	char Path[100] = "DATA";
	char PIfname[100] = "Pi.dat";
	char MUfname[100] = "Mu.dat";
	char Sfname[100] = "LTSigma.dat";
	char overmap[100] = "overMap.dat";
	char overbarmax[100] = "overBarMax.dat";
	char dataX[100] = "x.dat";
	char Nksizes[100] = "Nk.dat";
		
			
	errflag = 0;

/* READING PARAMETERS AND CHECKING THEM */
	
	while ((c = getopt(argc, argv, ":b:m:p:K:se:z:u:r:a:l:P:M:S:D:I:X:W:C:n:#:")) != -1) {
        switch(c) {
        case 'b': /* average overlap */
            BarOmega = strtod(optarg, pp);
            break;
        case 'm': /* maximum overlap */
            MaxOmega = strtod(optarg, pp);
            break;
		case 'p': /* # of dimensions */
            p = atoi(optarg);
            break;
        case 'K': /* # of components */
            K = atoi(optarg);
			break;
        case 's': /* spherical cluster */
            sph = 1;
			break;
        case 'e': /* maximum eccentricity */
            emax = strtod(optarg, pp);
			break;
        case 'z': /* smallest mixing proportion */
            PiLow = strtod(optarg, pp);
			break;
        case 'u': /* upper bound for Uniform */
            Ubound = strtod(optarg, pp);
			break;
        case 'r': /* # of resimulations */
			resN = atoi(optarg);
			break;
        case 'a': /*epsilon / accuracy */
            pars[0] = strtod(optarg, pp);
			pars[1] = strtod(optarg, pp);
			break;
        case 'l': /* max # of integration terms */
            lim = atoi(optarg);
			break;
        case 'P': /* Name of the file of mixing proportions */
			strcpy(PIfname, optarg);
			break;
        case 'M': /* Name of the file with means */
			strcpy(MUfname, optarg);
			break;
        case 'S': /* Name of the file with covariance matrices */
			strcpy(Sfname, optarg);
			break;
        case 'D': /* Name of the working directory */
			strcpy(Path, optarg);
			break;	    
        case 'I': /* Name of the file with true sample sizes */
			strcpy(Nksizes, optarg);
			break;
        case 'X': /* Name of the output file */
			strcpy(dataX, optarg);
			break;			
        case 'W': /* Name of the file with overlap map */
			strcpy(overmap, optarg);
			break;
        case 'C': /* Name of the file with overlap characteristics */
			strcpy(overbarmax, optarg);
			break;			
        case 'n': /* # of observations in simulated datasets */
            g = atoi(optarg);
			break;
        case '#': /* # of simulated mixtures */
            nf = atoi(optarg);
			break;

										       
        case ':': /* option with no operand */
            printf("Operand for option -%c is not specified correctly...\n", optopt);
			errflag = 1;
            break;
        case '?':
            printf("Unrecognized option: -%c\n", optopt);
			errflag = 2;
		}
	}

/* DELAING WITH WRONG INPUT */

	if (p < 1){
		printf("Operand for option -p is not specified correctly...\n");
		errflag = 1;
	}
	if (K < 1){
		printf("Operand for option -K is not specified correctly...\n");
		errflag = 1;
	}
	if ((sph != 0) & (sph != 1)){
		printf("Operand for option -s is not specified correctly...\n");
		errflag = 1;
	}
	if ((emax <= 0.0) | (emax >= 1.0)){
		printf("Operand for option -e is not specified correctly...\n");
		errflag = 1;
	}
	if ((PiLow <= 0.0) | (PiLow > 1.0)){
		printf("Operand for option -z is not specified correctly...\n");
		errflag = 1;
	}
	if (Ubound <= 0.0){
		printf("Operand for option -u is not specified correctly...\n");
		errflag = 1;
	}
	if (resN <= 0){
		printf("Operand for option -r is not specified correctly...\n");
		errflag = 1;
	}
	if (pars[0] <= 0.0){
		printf("Operand for option -a is not specified correctly...\n");
		errflag = 1;
	}
	if (lim <= 0){
		printf("Operand for option -l is not specified correctly...\n");
		errflag = 1;
	}
	if ((g < 0) & ((MaxOmega != 0.0) | (BarOmega != 0.0))){
		printf("Operand for option -n is not specified correctly...\n");
		errflag = 1;
	}
	if (nf < 1){
		printf("Operand for option -# is not specified correctly...\n");
		errflag = 1;
	}

	if (errflag != 0) exit(1);

	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, PIfname);
	strcpy(PIfname, buff);

	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, MUfname);
	strcpy(MUfname, buff);
	
	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, Sfname);
	strcpy(Sfname, buff);
	
	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, overmap);
	strcpy(overmap, buff);

	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, overbarmax);
	strcpy(overbarmax, buff);
	
	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, dataX);
	strcpy(dataX, buff);
	
	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, Nksizes);
	strcpy(Nksizes, buff);
			
	srand(time(NULL));

	MAKE_VECTOR(Pi, K);
	MAKE_MATRIX(Mu, K, p);
	MAKE_3ARRAY(S, K, p, p);

	MAKE_MATRIX(OmegaMap, K, K);
	MAKE_VECTOR(rcMax, 2);
			

/* RUNNING A CORRESPONDING FUNCTION */


	for (i=1; i <= nf; i++){	

		if (nf != 1){
			printf("\n\nMIXTURE MODEL #%i:\n\n", i);
		}
		
		if (BarOmega == 0.0){
			if (MaxOmega == 0.0){
				/* AVERAGE AND MAXIMUM OVERLAPS ARE NOT SPECIFIED */
				freadParameters(p, K, Pi, Mu, S, PIfname, MUfname, Sfname, i);
				ExactOverlap(p, K, Pi, Mu, S, pars, lim, OmegaMap, &BarOmega, &MaxOmega, rcMax);
				
				if ((BarOmega > MaxOmega) | (BarOmega*K*(K-1.0)/2 < MaxOmega) | (BarOmega < 0.0)){
					printf("Incorrect parameters. ");
					printf("Check the options -p and -K...\n");
					exit(2);
				}			
				
				printOverlap(K, OmegaMap, BarOmega, MaxOmega, rcMax);
				fprintOverlap(K, OmegaMap, BarOmega, MaxOmega, rcMax, i, overmap, overbarmax);
				BarOmega = 0.0;
				MaxOmega = 0.0;
				
				if (g > 0){
					double **Y;
					int *Nk;
					MAKE_MATRIX(Y, g, p);
					MAKE_VECTOR(Nk, K);
					
					/* dataset generation and printing parameters */
					genData(p, K, Pi, Mu, S, g, Y, Nk);
					fprintData(g, p, K, Y, Nk, i, dataX, Nksizes);
					printf("Dataset with cluster sizes Nk = ");
					for (k=0; k<K; k++){	
						printf("%i ",Nk[k]);
					}
					printf("has been generated...\n");
					
					FREE_MATRIX(Y);
					FREE_VECTOR(Nk);
				}
				
			} else {
				/* ONLY MAXIMUM OVERLAP IS SPECIFIED */				
				method = 1;
				OmegaClust(MaxOmega, method, p, K, PiLow, Ubound, emax, pars, lim, resN, sph, Pi, Mu, S, OmegaMap, &BarOmega, &MaxOmega, rcMax, &fail);
				
				if (fail != 0){
					printf("The desired overlap has not been met...\n");
					exit(2);
				}
				
				printf("The desired overlap has been met...\n");
				printOverlap(K, OmegaMap, BarOmega, MaxOmega, rcMax);
				printParameters(p, K, Pi, Mu, S);
				fprintOverlap(K, OmegaMap, BarOmega, MaxOmega, rcMax, i, overmap, overbarmax);
				fprintParameters(p, K, Pi, Mu, S, PIfname, MUfname, Sfname, i);
				BarOmega = 0.0;
				
				if (g != 0){
					double **Y;
					int *Nk;
					MAKE_MATRIX(Y, g, p);
					MAKE_VECTOR(Nk, K);
					
					/* dataset generation and printing parameters */
					genData(p, K, Pi, Mu, S, g, Y, Nk);
					fprintData(g, p, K, Y, Nk, i, dataX, Nksizes);
					printf("Dataset with cluster sizes Nk = ");
					for (k=0; k<K; k++){	
						printf("%i ",Nk[k]);
					}
					printf("has been generated...\n");
					
					FREE_MATRIX(Y);
					FREE_VECTOR(Nk);
				}
				
			}
		} else {
			if (MaxOmega == 0.0){
				/* ONLY AVERAGE OVERLAP IS SPECIFIED */				
				method = 0;
				OmegaClust(BarOmega, method, p, K, PiLow, Ubound, emax, pars, lim, resN, sph, Pi, Mu, S, OmegaMap, &BarOmega, &MaxOmega, rcMax, &fail);
				
				if (fail != 0){
					printf("The desired overlap has not been met...\n");
					exit(2);
				}
				
				printf("The desired overlap has been met...\n");
				printOverlap(K, OmegaMap, BarOmega, MaxOmega, rcMax);
				printParameters(p, K, Pi, Mu, S);
				fprintOverlap(K, OmegaMap, BarOmega, MaxOmega, rcMax, i, overmap, overbarmax);
				fprintParameters(p, K, Pi, Mu, S, PIfname, MUfname, Sfname, i);
				MaxOmega = 0.0;
				
				if (g != 0){
					double **Y;
					int *Nk;
					MAKE_MATRIX(Y, g, p);
					MAKE_VECTOR(Nk, K);
					
					/* dataset generation and printing parameters */
					genData(p, K, Pi, Mu, S, g, Y, Nk);
					fprintData(g, p, K, Y, Nk, i, dataX, Nksizes);
					printf("Dataset with cluster sizes Nk = ");
					for (k=0; k<K; k++){	
						printf("%i ",Nk[k]);
					}
					printf("has been generated...\n");					
					
					FREE_MATRIX(Y);
					FREE_VECTOR(Nk);
				}
				
			} else {
				/* AVERAGE AND MAXIMUM OVERLAPS ARE BOTH SPECIFIED */		
				OmegaBarOmegaMax(p, K, PiLow, Ubound, emax, pars, lim, resN, sph, Pi, Mu, S, OmegaMap, &BarOmega, &MaxOmega, rcMax, &fail);
				
				if (fail != 0){
					printf("The desired overlap has not been met...\n");
					exit(2);
				}
				
				printf("The desired overlap has been met...\n");
				printOverlap(K, OmegaMap, BarOmega, MaxOmega, rcMax);
				printParameters(p, K, Pi, Mu, S);
				fprintOverlap(K, OmegaMap, BarOmega, MaxOmega, rcMax, i, overmap, overbarmax);
				fprintParameters(p, K, Pi, Mu, S, PIfname, MUfname, Sfname, i);
				
				if (g != 0){
					double **Y;
					int *Nk;
					MAKE_MATRIX(Y, g, p);
					MAKE_VECTOR(Nk, K);
					
					/* dataset generation and printing parameters */
					genData(p, K, Pi, Mu, S, g, Y, Nk);
					fprintData(g, p, K, Y, Nk, i, dataX, Nksizes);
					printf("Dataset with cluster sizes Nk = ");
					for (k=0; k<K; k++){	
						printf("%i ",Nk[k]);
					}
					printf("has been generated...\n");
					
					FREE_MATRIX(Y);
					FREE_VECTOR(Nk);
				}
				
			}
		}

	}


	FREE_VECTOR(Pi);
	FREE_MATRIX(Mu);
	FREE_3ARRAY(S);

	FREE_MATRIX(OmegaMap);
	FREE_VECTOR(rcMax);

	FREE_VECTOR(pars);

	return 0;

}

