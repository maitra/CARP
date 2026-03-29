
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





/*
This is the main program that reads parameters and runs the following 3 steps:
	1. C-MixSim which generates mixtures with prespecified level of overlap	and datasets from these mixtures.
	2. Clustering algorithm provided by user.
	3. AdjRand which computes Adjusted Rand index for simulated datasets.
*/


#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<time.h>
#include<string.h>


int getopt(int argc, char *const argv[], const char *optstring);
extern char *optarg;
extern int optind, optopt;


int main(int argc, char *argv[])
{
	
	if (argc == 1){

		printf("Usage:\n> ./C-MixSim <set of parameters>\n\nParameters:\n");
		printf("-b: average overlap (no default value)\n");
		printf("-m: maximum overlap (no default value)\n");
		printf("-p: number of dimensions (2 by default)\n");
		printf("-K : number of mixing components (2 by default)\n");
		printf("-n: number of observations generated from every mixture (no default value)\n");
		printf("-#: number of simulated mixtures (1 by default)\n");
		printf("-0 : name of clustering program (no default name)\n");
		printf("-1 : name of partition analyzing program (“AdjRand” by default)\n");
		printf("-s: spherical covariance matrix structure (non-spherical by default if option is unspecified)\n");
		printf("-e: maximum eccentricity (0.90 by default)\n");
		printf("-z : smallest mixing proportion (equal mixing proportions 1/K by default)\n");
		printf("-u: upper bound for Uniform(0, upper-bound ) distribution from which mean vectors are generated\n");
		printf("-r : maximum number of resimulations (100 by default)\n");
		printf("-a: accuracy of estimation (1e-06 by default)\n");
		printf("-l : maximum number of integration terms (1e06 by default)\n");
		printf("-P : name of the file containing mixing proportions (“Pi.dat” by default)\n");
		printf("-M : name of the file containing mean vectors (“Mu.dat” by default)\n");
		printf("-S : name of the file containing covariance matrices in triangular form (“S.dat” by default)\n");
		printf("-D: name of working directory (“DATA” by default)\n");
		printf("-I : name of the file containing numbers of observations generated from every cluster (“Nk.dat” by default)\n");
		printf("-i : name of the file containing estimated classifications (“idEst.dat” by default)\n");
		printf("-X : name of the file containing simulated datasets (“x.dat” by default)\n");
		printf("-W : name of the file containing maps of pairwise overlaps (“overMap.dat” by default)\n");
		printf("-C : name of the file containing characteristics of simulated mixtures (“overBarMax.dat” by default)\n");
		printf("-R: name of the file containing index values (Adjusted Rand index and “AR.dat” by default)\n\n");

		exit(1);

	}	
	
	int c, errflag0, errflag1, code;

	char file1[500] = "./C-MixSim";
	char file2[500] = "./";
	char file3[500] = "./";
	char options1[500] = "";	
	char options2[500] = "";
	char options3[500] = "";

	errflag0 = 1;
	errflag1 = 1;
	
	/* DEALING WITH PARAMETERS */

	while ((c = getopt(argc, argv, ":b:m:p:K:se:z:u:r:a:l:P:M:S:D:I:i:X:W:C:R:n:#:0:1:")) != -1) {
		switch(c) {
        	case '0':
				strcat(file2, optarg);
				errflag0 = 0;
				break;
        	case '1':
				strcat(file3, optarg);
				errflag1 = 0;
				break;
			default:

				if ((c != '0') & (c != 'i') & (c != 'R')){
 					strcat(options1, " ");
					strcat(options1, argv[optind-1]);
				}

				if ((c == 'p') | (c == 'K') | (c == 'n') | (c == '#') | (c == 'D') | (c == 'X') | (c == 'i')){
					strcat(options2, " ");
					strcat(options2, argv[optind-1]);
				}
				
				if ((c == 'K') | (c == 'n') | (c == '#') | (c == 'D') | (c == 'I') | (c == 'i') | (c == 'R')){
					strcat(options3, " ");
					strcat(options3, argv[optind-1]);
				}				
				
				break;
		}
	}

/*	default analyzing program */

	if (errflag1 != 0){
		strcat(file3, "AdjRand");
	}

	strcat(file1, options1);
	strcat(file2, options2);
	strcat(file3, options3);

//	if (errflag0 != 0){
//		printf("ERROR: No clustering program has been detected [option -0]\n");
//		exit(1);
//	}
	
	
/*	STEP 1 - Mixture simulation	*/
	
	printf("\n\nSTEP 1: Running %s\n\n", file1);		
	code = system(file1);
	if (code != 0){
		printf("ERROR: Unable to fulfill STEP 1\n");
		exit(1);
	}

	if (errflag0 == 0){

/*		STEP 2 - Run of clustering algorithm	*/
	
		printf("\n\n\nSTEP 2: Running %s\n\n\n", file2);
		code = system(file2);
		if (code != 0){
			printf("ERROR: Unable to fulfill STEP 2\n");
			exit(1);
		}

/*		STEP 3 - Calculation of Adjusted Rand index	*/

		printf("\n\n\nSTEP 3: Running %s\n\n\n", file3);
		code = system(file3);
		if (code != 0){
			printf("ERROR: Unable to fulfill STEP 3\n");
			exit(1);
		}
	
	}
	
	return 0;
}

