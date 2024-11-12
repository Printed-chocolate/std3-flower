#ifndef STRAY_C
#define STRAY_C

#include <stdlib.h>
#include "nrutil.c"

#define PI 3.141592653589793

extern double dx,dy,dz;
extern int M,N,K,MM,NN,KK;

double stray(double ***m1,double ***m2,double ***m3,
	double ***hd1,double ***hd2,double ***hd3)
{
	int i,j,k;
	double stren;
	stren=0.0;
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				stren=stren+
					(hd1[i][j][k]*m1[i][j][k]+
					hd2[i][j][k]*m2[i][j][k]+
					hd3[i][j][k]*m3[i][j][k])*
					dx*dy*dz;
			}
		}
	}
	stren=-0.5*stren;
	return stren;

}	

#endif
	
	
