
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
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"
#define inf 1e+40;


/*  Given a hierarchical clustering through a sequence of agglomerations,
    derive the assignment for the top (lev-1)th level of the hierarchy. 
    
    Parameters:                                                    
    n:          number of observations                             
    ia,ib:      vectors of dimension N defining the agglomerations 
    lev:        number of clusters in partition.           
    iclass:     n-dimensional vector of cluster assignments       

    C code written by Ranjan Maitra, Baltimore (07/02/02) */ 

void hclass(int n, int *ia, int *ib, int lev, int *iclass)
{
  int i,j,k;
  for (i=0;i<n;i++) {
    iclass[i]=0;
  }
  j=lev-1;
  for (i=n-lev;i<(n-1);i++) {
    iclass[ib[i]]=j;
    for (k=(n-lev-1);k>=0;k--) {
      if (iclass[ia[k]]==j) {
	iclass[ib[k]]=j;
	}
    }
    j--;
  }
  j=0;
  iclass[ia[n-2]]=j;
  for (k=(n-lev-1);k>=0;k--) {
    if (iclass[ia[k]]==j) {
      iclass[ib[k]]=j;
    }
  }
  return;
}





/*  HIERARCHICAL CLUSTERING using (user-specified) criterion. 

    Parameters:                                               
    n,m               number of observations and coordinates  
    data(n,m)         input data matrix,                      
    iopt              clustering criterion to be used,        
    IA, IB, CRIT      history of agglomerations; dimension N-1 each       

    C code rewritten by Ranjan Maitra, Baltimore (07/02/02) from f77
    code of F. Murtagh, ESA/ESO/STECF, Garching (1986). Note that I use the 
    true distances, not squared Euclidean which FM used. */


int ioffset(int n, int i, int j)
{                               /*  Map (i,j)th element of upper half */
  return(j+i*n-(i+1)*(i+2)/2);  /*  diagonal symmetric matrix to vector */
} 

int imin(int aa, int bb)
{
  if (aa<bb) {
    return(aa);
  }
  else return(bb);
}

int imax(int aa, int bb)
{
  if (aa<bb) {
    return(bb);
  }
  else return(aa);
}

double fmin(double aa, double bb)
{
  if (aa<bb) {
    return(aa);
  }
  else return(bb);
}

double fmax(double aa, double bb)
{
  if (aa<bb) {
    return(bb);
  }
  else return(aa);
}

void hc(int n, int m, int iopt, double **data, int *ia, int *ib, double *crit)
{
  int i,i2,ind,ind1,ind2,ind3,ncl,j,j2,k,len,iflag,*nn,*flag;
  double dmin,r1,r2,x,xx,*disnn,*membr,*diss,*disso;

  ncl = n;
  len=n*(n-1)/2;
  
  MAKE_VECTOR(flag,n);
  MAKE_VECTOR(nn,n);
  MAKE_VECTOR(disnn,n);
  MAKE_VECTOR(membr,n);
  MAKE_VECTOR(disso,len);
  MAKE_VECTOR(diss,len);

  for (i=0; i<n; ++i) {
    nn[i]=0;
    disnn[i]=0;
    membr[i] = 1;
    flag[i] = 1;
  }

  for (i=0; i<len; ++i) {
    disso[i]=0;
    diss[i]=0;
  }

  for (i=0;i<(n-1);i++) {   /*  Construct dissimilarity matrix */
    for (j=(i+1);j<n;j++) {
      ind = ioffset(n, i, j);
      for (k=0;k<m;k++) { /* Computing Euclidean distance */
	diss[ind]+=(data[i][k]-data[j][k])*(data[i][k]-data[j][k]);
      }
      if (iopt == 1) {  /*for the case of the min. var. method where merging */
	diss[ind] /= 2.; /*criteria are defined in terms of variances */ 
      }
      else {
	diss[ind]=sqrt(diss[ind]); /* true Euclidean distance */
      }
      disso[ind] = diss[ind];
    }
  }
  
  for (i=0;i<(n-1);++i) { /*Carry out an agglomeration - first create list of NNs */
    int jm;
    dmin = inf;
    jm=i+1;
    for (j=(i+1);j<n;++j) {
      ind=ioffset(n,i,j);
      if (diss[ind]<dmin) {
	dmin = diss[ind];
	jm = j;
      }
    }
    nn[i]=jm;
    disnn[i]=dmin;
  }

  while (ncl>1) { /*     Next, determine least diss. using list of NNs */ 
    int jj=0,jm=0,im=0;
    dmin=inf;
    for (i= 0;i<(n-1);++i) {
      if (flag[i]==1) {
	if (disnn[i]<dmin) {
	  dmin=disnn[i];
	  im=i;
	  jm=nn[i];
	  }
      }
    }

    i2=imin(im,jm); /*  This allows an agglomeration to be carried out. */
    j2=imax(im,jm);
    ia[n-ncl]=i2;
    ib[n-ncl]=j2;
    crit[n-ncl]=dmin;
    
    flag[j2]=0;     /*  Update dissimilarities from new cluster. */
    dmin=inf;
    for (k=0;k<n;k++) {
      if ((flag[k]==1) && (k!=i2)) {
	x=membr[i2]+membr[j2]+membr[k];
	ind1=ioffset(n,imin(i2,k),imax(i2,k));
	ind2=ioffset(n,imin(j2,k),imax(j2,k));
	ind3=ioffset(n,i2,j2);
	xx=diss[ind3];
	
	if (iopt == 1) { /*  WARD'S MINIMUM VARIANCE METHOD - IOPT=1. */
	  diss[ind1] = (membr[i2] + membr[k]) * diss[ind1] + 
	    (membr[j2] + membr[k]) * diss[ind2] - membr[k] * xx;
	  diss[ind1] /= x;
	}
	
	if (iopt == 2) { /*  SINGLE LINK METHOD - IOPT=2. */ 
	  r1 = diss[ind1]; /* Computing MIN */
	  r2 = diss[ind2];
	  diss[ind1]=fmin(r1,r2);
	}
	
	if (iopt==3) {  /*  COMPLETE LINK METHOD - IOPT=3. */
	  r1=diss[ind1]; /* Computing MAX */
	  r2=diss[ind2];
	  diss[ind1] = fmax(r1,r2);
	}
	
	if (iopt==4) { /* AVG. LINK (OR GROUP AVERAGE) METHOD - IOPT=4. */
	  diss[ind1] = (membr[i2]*diss[ind1]+membr[j2]*diss[ind2])/
	    (membr[i2]+membr[j2]); 
	}
	
	if (iopt == 5) { /*  MCQUITTY'S METHOD - IOPT=5. */
	    diss[ind1] = diss[ind1]*0.5 + diss[ind2]*0.5;
	}
	
	if (iopt==6) { /*  MEDIAN (GOWER'S) METHOD - IOPT=6. */
	  diss[ind1]=diss[ind1]*0.5+diss[ind2]*0.5-xx*0.25;
	}
	
	if (iopt==7) { /*  CENTROID METHOD - IOPT=7. */
	  diss[ind1] = (membr[i2]*diss[ind1]+membr[j2]*diss[ind2] - 
			membr[i2]*membr[j2]*xx/(membr[i2]+membr[j2])
			)/(membr[i2] + membr[j2]);
	}
	if (i2<=k) {
	  if (diss[ind1]<dmin) { 
	    dmin=diss[ind1];
	    jj=k;
	  }
	}
      }
    }
    membr[i2]+=membr[j2];
    disnn[i2]=dmin;
    nn[i2] = jj;
    for (i=0;i<(n-1);++i) { /* Update list of NNs insofar as  required. */
      iflag=0;
      if (iflag==0) {
	if (flag[i]==0) {
	  iflag=1;
	}
	if (iflag==0) {
	  if ((nn[i]!=i2) && (nn[i]!=j2))  { 
	    iflag=1;
	  }
	  if (iflag==0) { /*  (Redetermine NN of i) */
	    int jj=0;
	    dmin = inf;
	    for (j=(i+1);j<n;++j) {
	      ind = ioffset(n,i,j);
	      if ((flag[j]==1) & (diss[ind]<=dmin)) {
		dmin = diss[ind];
		jj = j;
	      }
	    }
	    nn[i] = jj;
	    disnn[i] = dmin;
	  }
	}
      }
    }
    ncl--; /*  Repeat previous steps until n-1 agglomerations carried out. */
  }
  FREE_VECTOR(flag);
  FREE_VECTOR(nn);
  FREE_VECTOR(membr);
  FREE_VECTOR(disnn);
  FREE_VECTOR(disso);
  FREE_VECTOR(diss);
  return;
}



/* This code is a wrapper to the hc.c and hclass.c code. What it does is 
   return the class indicators for a hierarchically clustered tree, using some 
   user-specified criterion (out of a list of seven) and the number of classes 
   nclass */


void hclassify(int n,int m, double **x,int hcrit,int nclass,int *class)
{
  double *crit;
  int *ia,*ib;
  MAKE_VECTOR(ia,n);
  MAKE_VECTOR(ib,n);
  MAKE_VECTOR(crit,n);

  hc(n,m,hcrit,x,ia,ib,crit);
  FREE_VECTOR(crit);
  
  hclass(n,ia,ib,nclass,class);

  FREE_VECTOR(ia);
  FREE_VECTOR(ib);
  return;
}





int getopt(int argc, char *const argv[], const char *optstring);
extern char *optarg;
extern int optind, optopt;


/* This porgram obtains classifications for simulated datasets */


int main(int argc, char *argv[])
{
	
	int p, K, g, nf, hlink, k, i, j, c, errflag;
	int *id;
	double **X;
	
	FILE *fdata, *fout;		
	
	char buff[200] = "";
	char Path[100] = "DATA";
	char Xname[100] = "x.dat";
	char IDname[100] = "idEst.dat";
	
/* DEFAULT PARAMETERS */
	p = 2;
	K = 2;
	g = 0; /* sample size for generated datasets */
	nf = 1; /* # of simulated mixtures */
	
	hlink = 1; 	/*	1 - Ward's
					2 - Single
					3 - Complete
					4 - Average
					5 - McQuitty's
					6 - Median
					7 - Centroid
				*/
	
	errflag = 0;	
	
	/* DEALING WITH PARAMETERS */
	
	if (argc == 1){
		exit(1);
	}

	while ((c = getopt(argc, argv, ":p:K:n:D:i:X:#:")) != -1) {
        switch(c) {
		case 'p': /* # of dimensions */
            p = atoi(optarg);
            break;
        case 'K': /* # of components */
            K = atoi(optarg);
			break;
        case 'n': /* # of observations in datasets */
            g = atoi(optarg);
            g = fabs(g);
			break;
        case 'D': /* Name of the working directory */
			strcpy(Path, optarg);
			break;	    
        case 'i': /* Name for estimated classification */
			strcpy(IDname, optarg);
			break;
        case 'X': /* Name for the file with covariance matrices */
			strcpy(Xname, optarg);
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

	if (p < 1){
		printf("Operand for option -p is not specified correctly...\n");
		errflag = 1;
	}
	if (K < 1){
		printf("Operand for option -K is not specified correctly...\n");
		errflag = 1;
	}
	if (g < 0){
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
	strcat(buff, Xname);
	strcpy(Xname, buff);

	strcpy(buff, Path);
	strcat(buff, "/");
	strcat(buff, IDname);
	strcpy(IDname, buff);
	
	MAKE_MATRIX(X, g, p);
	MAKE_VECTOR(id, g);	
	
	
	if (!(fdata = fopen(Xname, "r"))){
		fprintf(stderr, "Error: could not open file %s...\n", Xname);
		exit(EXIT_FAILURE);
	}

	if (!(fout = fopen(IDname, "w"))){
		fprintf(stderr, "Error: could not open file %s...\n", IDname);
		exit(EXIT_FAILURE);
	}
	
	

	printf("OBTAINED CLASSIFICATIONS:\n");

	for (k=1; k<=nf; k++){

		if (nf != 1){
			printf("Dataset #%i:\n", k);
		}
		
	
		for (i=0; i<g; i++){
			for (j=0; j<p; j++){
				if(!fscanf(fdata, "%lf ", &X[i][j])) {
					fprintf(stderr, "Error: invalid format in %s...\n", Xname);
					exit(EXIT_FAILURE);
				}
			}
		}
	
		
		/* Run hierarchical clustering */
		
		hclassify(g, p, X, hlink, K, id);
	
		for (i=0; i<g; i++) {
			printf("%i ", id[i]);
		}
		printf("\n");
	

		for (i=0; i<g; i++) {
			fprintf(fout, "%i ", id[i]);
		}
	
	}

	fclose(fout);
	fclose(fdata);
	
	FREE_MATRIX(X);
	FREE_VECTOR(id);	

	return 0;

}

