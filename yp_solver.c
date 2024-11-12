#ifndef YP_SOLVER_C
#define YP_SOLVER_C
#include <math.h>
#include "nrutil.c"
#include "realft.c"
#define PI 3.141592653589793

void yp_solver(double ***u,double ***u1,double ***F,
	unsigned long  M, unsigned long N, unsigned long K,
	double Ddtx,double Ddty,double Ddtz,double dt)
{
	void realft(double *data,unsigned long n,int isign);

	double ***fp,***up,***up1;
	double tmp1,tmp2,tmp3,tmp4;	
	int i,j,k,j1;
	fp=f3tensor(1,M,1,K,1,N);
	up=f3tensor(1,M,1,K,1,N);
	up1=f3tensor(1,M,1,K,1,N);
	for (i=1;i<=M;i++){
		for (k=1;k<=K;k++){
			for (j=1;j<=N;j++){
				fp[i][k][j]=F[i][j][k];
				up[i][k][j]=u[i][j][k];
			}
		} 
	}
	for (i=1;i<=M;i++){
		for (k=1;k<=K;k++){
			realft(up[i][k],N,1);
			realft(fp[i][k],N,1);
		}
	}
	
	for (i=1;i<=M;i++){
		for (k=1;k<=K;k++){
			for (j1=1;j1<=0.5*N;j1++){	
				j=2*j1-1;
				if(i>1)tmp1=up[i-1][k][j];
				else tmp1=up[1][k][j];
				if(i<M)tmp2=up[i+1][k][j];
				else tmp2=up[M][k][j];
				if(k>1)tmp3=up[i][k-1][j];
				else tmp3=up[i][1][j];
				if(k<K)tmp4=up[i][k+1][j];
				else tmp4=up[i][K][j];	
				up1[i][k][j]=(up[i][k][j]*
					(1.0-2.0*(Ddtx+Ddtz))+
					Ddtx*(tmp1+tmp2)+Ddtz*(tmp3+tmp4)
					+dt/3.0*fp[i][k][j])/
					(1.0-2.0*Ddty*(cos(2.0*PI*(j1-1.0)/N)-1.0));
				j=2*j1;
				if(i>1)tmp1=up[i-1][k][j];
				else tmp1=up[1][k][j];
				if(i<M)tmp2=up[i+1][k][j];
				else tmp2=up[M][k][j];
				if(k>1)tmp3=up[i][k-1][j];
				else tmp3=up[i][1][j];
				if(k<K)tmp4=up[i][k+1][j];
				else tmp4=up[i][K][j];	
				up1[i][k][j]=(up[i][k][j]*
					(1.0-2.0*(Ddtx+Ddtz))+
					Ddtx*(tmp1+tmp2)+Ddtz*(tmp3+tmp4)
					+dt/3.0*fp[i][k][j])/
					(1.0-2.0*Ddty*(cos(2.0*PI*(j1-1.0)/N)-1.0));
			}
			up1[i][k][2]=up1[i][k][2]/(1.0+4.0*Ddty);
		}
	}
	for (i=1;i<=M;i++){
		for (k=1;k<=K;k++){
			realft(up1[i][k],N,-1);
				for (j=1;j<=N;j++){
					up1[i][k][j]=2.0*up1[i][k][j]/N;
				}
		}
	}
	for (i=1;i<=M;i++){
		for(k=1;k<=K;k++){
			for (j=1;j<=N;j++){
				u1[i][j][k]=up1[i][k][j];
		   }
		}
	}
	free_f3tensor(fp,1,M,1,K,1,N);
	free_f3tensor(up,1,M,1,K,1,N);
	free_f3tensor(up1,1,M,1,K,1,N);

}
#undef PI
#endif
	
