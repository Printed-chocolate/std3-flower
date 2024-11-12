#ifndef XP_SOLVER_C
#define XP_SOLVER_C
#include <math.h>
#include "nrutil.c"
#include "realft.c"
#define PI 3.141592653589793

void xp_solver(double ***u,double ***u1,double ***F,
	unsigned long  M, unsigned long N, unsigned long K,
	double Ddtx,double Ddty,double Ddtz,double dt)
{
	void realft(double *data,unsigned long n,int isign);

	double ***fp,***up,***up1;
	double tmp1,tmp2,tmp3,tmp4;	
	int i,j,k,i1;
	fp=f3tensor(1,N,1,K,1,M);
	up=f3tensor(1,N,1,K,1,M);
	up1=f3tensor(1,N,1,K,1,M);
	for (j=1;j<=N;j++){
		for (k=1;k<=K;k++){
			for (i=1;i<=M;i++){
				fp[j][k][i]=F[i][j][k];
				up[j][k][i]=u[i][j][k];
			}
		}
	}
	for (j=1;j<=N;j++){
		for (k=1;k<=K;k++){
			realft(up[j][k],M,1);
			realft(fp[j][k],M,1);
		}
	}
	
	for (j=1;j<=N;j++){
		for (k=1;k<=K;k++){
			for (i1=1;i1<=0.5*M;i1++){
				i=2*i1-1;
				if(j>1)tmp1=up[j-1][k][i];
				else tmp1=up[1][k][i];
				if(j<N)tmp2=up[j+1][k][i];
				else tmp2=up[N][k][i];
				if(k>1)tmp3=up[j][k-1][i];
				else tmp3=up[j][1][i];
				if(k<K)tmp4=up[j][k+1][i];
				else tmp4=up[j][K][i];	
				up1[j][k][i]=(up[j][k][i]*
					(1.0-2.0*(Ddty+Ddtz))+
					Ddty*(tmp1+tmp2)+Ddtz*(tmp3+tmp4)
					+dt/3.0*fp[j][k][i])/
					(1.0-2.0*Ddtx*(cos(PI*(i1-1.0)/M)-1.0));
				i=2*i1;
				if(j>1)tmp1=up[j-1][k][i];
				else tmp1=up[1][k][i];
				if(j<N)tmp2=up[j+1][k][i];
				else tmp2=up[N][k][i];
				if(k>1)tmp3=up[j][k-1][i];
				else tmp3=up[j][1][i];
				if(k<K)tmp4=up[j][k+1][i];
				else tmp4=up[j][K][i];	
				up1[j][k][i]=(up[j][k][i]*
					(1.0-2.0*(Ddty+Ddtz))+
					Ddty*(tmp1+tmp2)+Ddtz*(tmp3+tmp4)
					+dt/3.0*fp[j][k][i])/
					(1.0-2.0*Ddtx*(cos(PI*(i1-1.0)/M)-1.0));
			}
			up1[j][k][2]=up1[j][k][2]/(1.0+4.0*Ddtx);
		}
	}
	for (j=1;j<=N;j++){
		for (k=1;k<=K;k++){
			realft(up1[j][k],M,-1);
				for (i=1;i<=M;i++){
					up1[j][k][i]=2.0*up1[j][k][i]/M;
				}
		}
	}
	for (j=1;j<=N;j++){
		for(k=1;k<=K;k++){
			for (i=1;i<=M;i++){
				u1[i][j][k]=up1[j][k][i];
		   }
		}
	}
	free_f3tensor(fp,1,M,1,K,1,N);
	free_f3tensor(up,1,M,1,K,1,N);
	free_f3tensor(up1,1,M,1,K,1,N);

}
#undef PI
#endif
	
