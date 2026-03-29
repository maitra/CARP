
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

/* fprintOverlap : writes overlaps to files
 * Parameters:
 * 		K - number of components
 * 		OmegaMap  - map of misclassification probabilities
 *  	BarOmega - average overlap
 * 		MaxOmega - maximum overlap
 * 		rcMax - pair of components producing maximum overlap
 * 		fileN - number of the dataset
 * 		overmap - name of file with the map of misclassification probabilities
 * 		overbarmax - name of file with overlap characteristics
 */

void fprintOverlap(int K, double **OmegaMap, double BarOmega, double MaxOmega, int *rcMax, int fileN, char *overmap, char *overbarmax){

	int i, j;
	
	FILE *fout;

	if (fileN == 1){
		if (!(fout = fopen(overmap, "w"))){
			fprintf(stderr, "Error: could not open file %s...\n", overmap);
			exit(EXIT_FAILURE);
		}
	} else {
		if (!(fout = fopen(overmap, "a"))){
			fprintf(stderr, "Error: could not open file %s...\n", overmap);
			exit(EXIT_FAILURE);
		}
	}

	for (i=0; i<K; i++){
		for (j=0; j<K; j++){
			fprintf(fout, "%lf ", OmegaMap[i][j]);
		}
		fprintf(fout, "\n");
	}

	fclose(fout);


	if (fileN == 1){
		if (!(fout = fopen(overbarmax, "w"))){
			fprintf(stderr, "Error: could not open file %s...\n", overbarmax);
			exit(EXIT_FAILURE);
		}
	} else {
		if (!(fout = fopen(overbarmax, "a"))){
			fprintf(stderr, "Error: could not open file %s...\n", overbarmax);
			exit(EXIT_FAILURE);
		}
	}

	fprintf(fout, "%lf %lf %i %i\n", BarOmega, MaxOmega, rcMax[0], rcMax[1]);

	fclose(fout);

}


/* printOverlap : prints overlaps
 * Parameters:
 * 		K - number of components
 * 		OmegaMap - map of misclassification probabilities
 * 		BarOmega - average overlap
 * 		MaxOmega - maximum overlap
 * 		rcMax - pair of components producing maximum overlap
 */

void printOverlap(int K, double **OmegaMap, double BarOmega, double MaxOmega, int *rcMax){

	int i, j;

	printf("Map of misclassification probabilities:\n");
	for (i=0; i<K; i++){		
		for (j=0; j<K; j++){
			printf("[%i][%i]:%lf ", i, j, OmegaMap[i][j]);
		}
		printf("\n");
	}

	printf("Average Overlap: %lf\n", BarOmega);
	printf("Maximum Overlap: %lf", MaxOmega);
	printf(" (components: %i and %i)\n", rcMax[0], rcMax[1]);

}


/* fprintParameters : writes parameters to files
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		Pi - vector of mixing proportions
 * 		Mu - set of mean vectors
 * 		S - set of covariance matrices
 * 		PIfname - name of file with mixing proportions
 * 		MUfname - name of file with mean vectors
 * 		Sfname - name of file with covariance matrices
 * 		fileN - number of the dataset
 */

void fprintParameters(int p, int K, double *Pi, double **Mu, double ***S, char *PIfname, char *MUfname, char *Sfname, int fileN){

	int i, j, k;

	FILE *fout;

	if (fileN == 1){
		if (!(fout = fopen(PIfname, "w"))){
			fprintf(stderr, "Error: could not open file %s...\n", PIfname);
			exit(EXIT_FAILURE);
		}
	} else {
		if (!(fout = fopen(PIfname, "a"))){
			fprintf(stderr, "Error: could not open file %s...\n", PIfname);
			exit(EXIT_FAILURE);
		}
	}

	for (k=0; k<K; k++){		
		fprintf(fout, "%e ", Pi[k]);
	}
	fprintf(fout, "\n");	
	
	fclose(fout);


	if (fileN == 1){
		if (!(fout = fopen(MUfname, "w"))){
			fprintf(stderr, "Error: could not open file %s...\n", MUfname);
			exit(EXIT_FAILURE);
		}
	} else {
		if (!(fout = fopen(MUfname, "a"))){
			fprintf(stderr, "Error: could not open file %s...\n", MUfname);
			exit(EXIT_FAILURE);
		}
	}

	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
			fprintf(fout, "%e ", Mu[k][i]);
		}
		fprintf(fout, "\n");
	}

	fclose(fout);


	if (fileN == 1){
		if (!(fout = fopen(Sfname, "w"))){
			fprintf(stderr, "Error: could not open file %s...\n", Sfname);
			exit(EXIT_FAILURE);
		}
	} else {
		if (!(fout = fopen(Sfname, "a"))){
			fprintf(stderr, "Error: could not open file %s...\n", Sfname);
			exit(EXIT_FAILURE);
		}
	}

	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
			for (j=0; j<(i+1); j++){
				fprintf(fout, "%e ", S[k][i][j]);
			}
		}
		fprintf(fout, "\n");
	}

	fclose(fout);

}


/* prints the parameters
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		Pi - vector of mixing proportions
 * 		Mu - set of mean vectors
 * 		S - set of covariance matrices
 */

void printParameters(int p, int K, double *Pi, double **Mu, double ***S){

	int i, j, k;

	printf("\nMixture parameters:\n");
	printf("Pi:\n");
	for (k=0; k<K; k++){		
		printf("%lf ",Pi[k]);
	}
	printf("\n");
	printf("Mu:\n");
	for (k=0; k<K; k++){
		printf("[%i]: ",k);
		for (i=0; i<p; i++){
			printf("%lf ", Mu[k][i]);
		}
		printf("\n");
	}
	printf("Sigma:\n");
	for (k=0; k<K; k++){
		printf("[%i]:\n",k);
		for (i=0; i<p; i++){
			printf("  ");
			for (j=0; j<p; j++){
				printf("%lf ", S[k][i][j]);
			}
			printf("\n");
		}
	}

}


/* freadParameters: reads the parameters from files
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		Pi - vector of mixing proportions
 * 		Mu - set of mean vectors
 * 		S - set of covariance matrices
 * 		PIfname - name of file with mixing proportions
 * 		MUfname - name of file with mean vectors
 * 		Sfname - name of file with covariance matrices
 * 		fileN - number of the dataset
 */
 
void freadParameters(int p, int K, double *Pi, double **Mu, double ***S, char *PIfname, char *MUfname, char *Sfname, int fileN){

	int i, j, k, r;

	FILE *finp;

	if (!(finp = fopen(PIfname, "r"))){
		fprintf(stderr, "Error: could not open file %s...\n", PIfname);

		exit(EXIT_FAILURE);
        }

	for (r=1; r<=fileN; r++){
		for (k=0; k<K; k++){
			if(!fscanf(finp, "%lf\n", &Pi[k])){
				fprintf(stderr, "Error: invalid format in %s...\n", PIfname);
				exit(EXIT_FAILURE);
			}
		}
	}



	if (!(finp = fopen(MUfname, "r"))){
		fprintf(stderr, "Error: could not open file %s...\n", MUfname);

		exit(EXIT_FAILURE);
        }

	for (r=1; r<=fileN; r++){
		for (k=0; k<K; k++){
			for (i=0; i<p; i++){
				if(!fscanf(finp, "%lf ", &Mu[k][i])) {
					fprintf(stderr, "Error: invalid format in %s...\n", MUfname);
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	fclose(finp);



	if (!(finp = fopen(Sfname, "r"))){
		fprintf(stderr, "Error: could not open file %s...\n", Sfname);
		exit(EXIT_FAILURE);
        }

	for (r=1; r<=fileN; r++){
		for (k=0; k<K; k++){
			for (i=0; i<p; i++){
				for (j=0; j<(i+1); j++){
					if(!fscanf(finp, "%lf ", &S[k][i][j])) {
						fprintf(stderr, "Error: invalid format in %s...\n", Sfname);
						exit(EXIT_FAILURE);
					}
					if (i != j) S[k][j][i] = S[k][i][j];
				}
			}
		}
	}

	fclose(finp);

}


/* fprintData : writes simulated datasets and cluster sizes to files
 * Parameters:
 * 		n - number of observations
 * 		p - number of dimensions
 * 		K - number of components
 * 		x - dataset
 * 		Nk - cluster sizes
 * 		fileN - number of the dataset
 * 		dataX - name of the output data file
 * 		Nksizes - name of the output classification file
 */

void fprintData(int n, int p, int K, double **x, int *Nk, int fileN, char *dataX, char *Nksizes){

	int i, j, k;
	
	FILE *fout;

	if (fileN == 1){
		if (!(fout = fopen(dataX, "w"))){
			fprintf(stderr, "Error: could not open file %s...\n", dataX);
			exit(EXIT_FAILURE);
		}
	} else {
		if (!(fout = fopen(dataX, "a"))){
			fprintf(stderr, "Error: could not open file %s...\n", dataX);
			exit(EXIT_FAILURE);
		}
	}

	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			fprintf(fout, "%lf ", x[i][j]);
		}
		fprintf(fout, "\n");
	}

	fclose(fout);

	if (fileN == 1){
		if (!(fout = fopen(Nksizes, "w"))){
			fprintf(stderr, "Error: could not open file %s...\n", Nksizes);
			exit(EXIT_FAILURE);
		}
	} else {
		if (!(fout = fopen(Nksizes, "a"))){
			fprintf(stderr, "Error: could not open file %s...\n", Nksizes);
			exit(EXIT_FAILURE);
		}
	}

	for (k=0; k<K; k++){
		fprintf(fout, "%i ", Nk[k]);
       	}
	fprintf(fout, "\n");

	fclose(fout);


}


