#ifndef ZSOLVER_C
#define ZSOLVER_C
#include "nrutil.c"
#include "tridag.c"
#define PI 3.141592653589793

void zsolver(double ***u,double ***u1,double ***f,
	unsigned long  M, unsigned long N, unsigned long K,
	double Ddtx,double Ddty,double Ddtz,double dt)
{
	void tridag(double a, double b, double c,
	double b1, double *r, double *u,
	unsigned long n);
	
	double *r,*utmp;
	double tmp1,tmp2,tmp3,tmp4;
	int i,j,k;
	r=vector(1,K);
	utmp=vector(1,K);
	
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){	
				if(i>1)tmp1=u[i-1][j][k];
				else tmp1=u[1][j][k];
				if(i<M)tmp2=u[i+1][j][k];
				else tmp2=u[M][j][k];
				if(j>1)tmp3=u[i][j-1][k];
				else tmp3=u[i][1][k];
				if(j<N)tmp4=u[i][j+1][k];
				else tmp4=u[i][N][k];	
				r[k]=u[i][j][k]*(1.0-2.0*(Ddtx+Ddty))+
					Ddtx*(tmp1+tmp2)+Ddty*(tmp3+tmp4)+
					dt/3.0*f[i][j][k];
			}
			tridag(-Ddtz,1.0+2.0*Ddtz,-Ddtz,1.0+Ddtz,r,utmp,K);
			for (k=1;k<=K;k++) u1[i][j][k]=utmp[k];
	
		}
	}
   free_vector(r,1,K);
   free_vector(utmp,1,K);
}
#undef PI
#endif