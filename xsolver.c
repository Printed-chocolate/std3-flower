#ifndef XSOLVER_C
#define XSOLVER_C
#include <math.h>
#include "nrutil.c"
#include "tridag.c"
#define PI 3.141592653589793
void xsolver(double ***u,double ***u1,double ***f,
	unsigned long  M, unsigned long N, unsigned long K,
	double Ddtx,double Ddty,double Ddtz,double dt)
{
	void tridag(double a, double b, double c,
	double b1, double *r, double *u,
	unsigned long n);
	
	double *r,*utmp;
	double tmp1,tmp2,tmp3,tmp4;
	int i,j,k;
	r=vector(1,M);
	utmp=vector(1,M);
	
	for (j=1;j<=N;j++){
		for (k=1;k<=K;k++){
			for (i=1;i<=M;i++){	
				if(j>1)tmp1=u[i][j-1][k];
				else tmp1=u[i][1][k];
				if(j<N)tmp2=u[i][j+1][k];
				else tmp2=u[i][N][k];
				if(k>1)tmp3=u[i][j][k-1];
				else tmp3=u[i][j][1];
				if(k<K)tmp4=u[i][j][k+1];
				else tmp4=u[i][j][K];	
				r[i]=u[i][j][k]*(1.0-2.0*(Ddty+Ddtz))+
					Ddty*(tmp1+tmp2)+Ddtz*(tmp3+tmp4)+
					dt/3.0*f[i][j][k];
			}
			tridag(-Ddtx,1.0+2.0*Ddtx,-Ddtx,1.0+Ddtx,r,utmp,M);
			for (i=1;i<=M;i++) u1[i][j][k]=utmp[i];
	
		}
	}
   free_vector(r,1,M);
   free_vector(utmp,1,M);
}
#undef PI
#endif
	
