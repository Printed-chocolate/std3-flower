#ifndef YSOLVER_C
#define YSOLVER_C
#include "nrutil.c"
#include "tridag.c"
#define PI 3.141592653589793

void ysolver(double ***u,double ***u1,double ***f,
	unsigned long  M, unsigned long N, unsigned long K,
	double Ddtx,double Ddty,double Ddtz,double dt)
{
	void tridag(double a, double b, double c,
	double b1, double *r, double *u,
	unsigned long n);
	
	double *r,*utmp;
	double tmp1,tmp2,tmp3,tmp4;
	int i,j,k;
	r=vector(1,N);
	utmp=vector(1,N);
	
	for (i=1;i<=M;i++){
		for (k=1;k<=K;k++){
			for (j=1;j<=N;j++){	
				if(i>1)tmp1=u[i-1][j][k];
				else tmp1=u[1][j][k];
				if(i<M)tmp2=u[i+1][j][k];
				else tmp2=u[M][j][k];
				if(k>1)tmp3=u[i][j][k-1];
				else tmp3=u[i][j][1];
				if(k<K)tmp4=u[i][j][k+1];
				else tmp4=u[i][j][K];	
				r[j]=u[i][j][k]*(1.0-2.0*(Ddtx+Ddtz))+
					Ddtx*(tmp1+tmp2)+Ddtz*(tmp3+tmp4)+
					dt/3.0*f[i][j][k];
			}
			tridag(-Ddty,1.0+2.0*Ddty,-Ddty,1.0+Ddty,r,utmp,N);
			for (j=1;j<=N;j++) u1[i][j][k]=utmp[j];
	
		}
	}
   free_vector(r,1,N);
   free_vector(utmp,1,N);
}
#undef PI
#endif
	
