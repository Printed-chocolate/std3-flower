#ifndef Z_SOLVER_C
#define Z_SOLVER_C
#include <math.h>
#include "nrutil.c"
#include "cosft2.c"
#define PI 3.141592653589793

void z_solver(double ***u,double ***u1,double ***F,
	unsigned long  M, unsigned long N, unsigned long K,
	double Ddtx,double Ddty,double Ddtz,double dt)
{
	double ***fp,***up,***up1;
	double tmp1,tmp2,tmp3,tmp4;
	int i,j,k;
	fp=f3tensor(1,M,1,N,1,K);
	up=f3tensor(1,M,1,N,1,K);
	up1=f3tensor(1,M,1,N,1,K);
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				fp[i][j][k]=F[i][j][k];
				up[i][j][k]=u[i][j][k];
			}
		}
	}
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			cosft2(up[i][j],K,1);
			cosft2(fp[i][j],K,1);
		}
	}
	
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				if(i>1)tmp1=up[i-1][j][k];
				else tmp1=up[1][j][k];
				if(i<M)tmp2=up[i+1][j][k];
				else tmp2=up[M][j][k];
				if(j>1)tmp3=up[i][j-1][k];
				else tmp3=up[i][1][k];
				if(j<N)tmp4=up[i][j+1][k];
				else tmp4=up[i][N][k];	
				up1[i][j][k]=(up[i][j][k]*
					(1.0-2.0*(Ddtx+Ddty))+
					Ddtx*(tmp1+tmp2)+Ddty*(tmp3+tmp4)
					+dt/3.0*fp[i][j][k])/
					(1.0-2.0*Ddtz*(cos(PI*(k-1.0)/K)-1.0));
			}
		}
	}
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			cosft2(up1[i][j],K,-1);
				for (k=1;k<=K;k++){
					up1[i][j][k]=2.0*up1[i][j][k]/K;
				}
		}
	}
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for(k=1;k<=K;k++){	
				u1[i][j][k]=up1[i][j][k];
			}
		}
	}
	free_f3tensor(fp,1,M,1,N,1,K);
	free_f3tensor(up,1,M,1,N,1,K);
	free_f3tensor(up1,1,M,1,N,1,K);

}
#undef PI
#endif
	
