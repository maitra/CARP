
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
#include<string.h>
#include<math.h>
#include "array.h"
#define inf 1e+40;

/* This procedure computes Adjusted Rand index
 * 
 *  Parameters:
 *  N - number of observations
 * 	TRUK - true number of clusters
 * 	PREDK - estimated number of clusters
 * 	trcl - true classification vector
 * 	prcl - estimated classification vector
 * 	Rand - Rand index
 * 	adjRand - Adjusted Rand index
 * 	Eindex - E index
 */

void RRand(int N, int TRUK, int PREDK,int *trcl, int *prcl,
	   double *Rand,double *adjRand,double *Eindex)
{
  int i,j,n[TRUK][PREDK];
  double sumtr[TRUK],sumpr[PREDK],sumprsq,sumtrsq,sumsq,discordant,
    sumtrprsq;
  double term1, term2, term3;
  double nij2sum, nidot2sum, ndotj2sum, Wallace;

  for (i=0;i<TRUK;i++)    {
    for (j=0;j<PREDK;j++)     {
      n[i][j]=0;
    }
  }
  
  for (i=0;i<N;i++) {
    n[trcl[i]][prcl[i]]+=1;
  }

  sumtrsq=0.;
  for (i=0;i<TRUK;i++) {
    sumtr[i]=0.;
    for (j=0;j<PREDK;j++) {
      sumtr[i]+=n[i][j];    }
    sumtrsq+=sumtr[i]*sumtr[i];
  }
  
  sumprsq=0.;
  for (j=0;j<PREDK;j++) {
    sumpr[j]=0.;
    for (i=0;i<TRUK;i++) {
      sumpr[j]+=(double)n[i][j];    }
    sumprsq+=sumpr[j]*sumpr[j];
  }

  sumtrprsq=0.;
  for (i=0;i<TRUK;i++) {
    for (j=0;j<PREDK;j++) {
      sumtrprsq+=sumtr[i]*sumtr[i]*sumpr[j]*sumpr[j];
    }
  }

  (*Eindex)=sumtrprsq/(N*((double)N-1) + N*(double)N/(N-1)) - (sumprsq + sumtrsq)/(N-1);
  (*Eindex)*=2.;
  (*Eindex)/=N*((double)N-1);
  
  sumsq=0.;
  for (i=0;i<TRUK;i++)    {
    for (j=0;j<PREDK;j++) {
      sumsq+=(double)n[i][j]*n[i][j];
    }
  }

  nij2sum=0.;
  for (i=0;i<TRUK;i++) {
    for (j=0;j<PREDK;j++) {
      nij2sum+=(double)n[i][j]*(n[i][j]-1)/2.0;
    }
  }

  nidot2sum=0.;
  for (i=0;i<TRUK;i++) {
    nidot2sum+=(double)sumtr[i]*(sumtr[i]-1)/2.0;
  }

  ndotj2sum=0.;
  for (i=0;i<PREDK;i++) {
    ndotj2sum+=(double)sumpr[i]*(sumpr[i]-1)/2.0;
  }

  Wallace=nij2sum/nidot2sum;
  discordant = 0.5*(sumtrsq + sumprsq) - sumsq ;

  (*Rand)=1.0-discordant/((double)N*((double)N-1.)/2.);

  term3 = nidot2sum * ndotj2sum / ((double)N*((double)N-1.)/2.);

  term1 = nij2sum - term3;

  term2 = (nidot2sum + ndotj2sum)/2 - term3;

  (*adjRand)= term1/term2;

}


/* This program computes Adjusted Rand index for simulated datasets
 */


int getopt(int argc, char *const argv[], const char *optstring);
extern char *optarg;
extern int optind, optopt;


int main(int argc, char *argv[])
{

	int K, g, nf, nk, k, i, j, c, start, errflag;
	int *idtrue, *idest;
	double Rand, adjRand, Eindex;
	
	FILE *fidtrue, *fidest, *fout;
	
	char buff[200] = "";
	char Path[100] = "DATA";
	char IDtrue[100] = "Nk.dat";
	char IDest[100] = "idEst.dat";
	char ARname[100] = "AR.dat";

/* DEFAULT PARAMETERS */
	K = 2;
	g = 0; /* sample size for generated datasets */
	nf = 1; /* # of simulated mixtures */
	
	errflag = 0;


	if (argc == 1){
		exit(1);
	}

	while ((c = getopt(argc, argv, ":K:n:D:I:i:R:#:")) != -1) {
        switch(c) {
        case 'K': /* # of components */
            K = atoi(optarg);
			break;
        case 'n': /* # of observations in simulated datasets */
            g = atoi(optarg);
            g = fabs(g);
			break;
        case 'D': /* Name working directory */
			strcpy(Path, optarg);
			break;	   			
        case 'I': /* true classification */
			strcpy(IDtrue, optarg);
			break;
        case 'i': /* estimated classification */
			strcpy(IDest, optarg);
			break;
        case 'R': /* AR index */
			strcpy(ARname, optarg);
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


	if (K < 1){
		printf("Operand for option -K is not specified correctly...\n");
		errflag = 1;
	}
	if (g <= 0){
		printf("Operand for option -g is not specified correctly...\n");
		errflag = 1;
	}
	if (nf < 1){
		printf("Operand for option -# is not specified correctly...\n");
		errflag = 1;
	}

	if (errflag != 0) exit(1);

	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, IDtrue);
	strcpy(IDtrue, buff);

	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, IDest);
	strcpy(IDest, buff);
	
	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, ARname);
	strcpy(ARname, buff);	
	
	MAKE_VECTOR(idtrue, g);
	MAKE_VECTOR(idest, g);	
	
	fidtrue = fopen(IDtrue,"r");
	fidest = fopen(IDest,"r");
	fout = fopen(ARname,"w");

	printf("ADJUSTED RAND VALUES:\n");

	for (k=1; k<=nf; k++){

		if (nf != 1){
			printf("Dataset #%i: ", k);
		}
		
		for (i=0; i<g; i++){
			fscanf(fidest, "%i ", &idest[i]);
		}

		start = 0;
		for (j=0; j<K; j++){
			fscanf(fidtrue, "%i ", &nk);
			for (i=start; i<(nk+start); i++){
				idtrue[i] = j;
			}
			start = start + nk;
		}		
	

		RRand(g, K, K, idtrue, idest, &Rand, &adjRand, &Eindex);
	
		printf("%lf\n", adjRand);
	
		fprintf(fout, "%lf ", adjRand);
	
	}

	fclose(fout);
	fclose(fidtrue);
	fclose(fidest);	

	FREE_VECTOR(idtrue);
	FREE_VECTOR(idest);	

	return 0;

}


