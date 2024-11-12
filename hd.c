#ifndef HD_C
#define HD_C

#include <stdlib.h>
#include "nrutil.c"
#include "conv3.c"
#include "kernel.c"

extern double dx,dy,dz,mspi;
extern int M,N,K,MM,NN,KK;

void hd(double ***m1,double ***m2,double ***m3,
	double ***hd1,double ***hd2,double ***hd3)
{
	void conv3(double ***data1,double **speq1,double ***data2,
		unsigned long n1,unsigned long n2,unsigned long n3);
	void getkernel(double ***kxx,double ***kyy,double ***kzz,
		double ***kxy,double ***kyz,double ***kzx,
		double **spxx,double **spyy,double **spzz,
		double **spxy,double **spyz,double **spzx);
	
	double ***mm1,***mm2,***mm3;
	double ***kxx,***kyy,***kzz,***kxy,***kyz,***kzx,
			***kyx,***kzy,***kxz;
	double **spxx,**spyy,**spzz,**spxy,**spyz,**spzx,
			**spyx,**spzy,**spxz;
	int i,j,k;
	
	mm1=f3tensor(1,MM,1,NN,1,KK);
	mm2=f3tensor(1,MM,1,NN,1,KK);
	mm3=f3tensor(1,MM,1,NN,1,KK);
	kxx=f3tensor(1,MM,1,NN,1,KK);
	kyy=f3tensor(1,MM,1,NN,1,KK);
	kzz=f3tensor(1,MM,1,NN,1,KK);
	kxy=f3tensor(1,MM,1,NN,1,KK);
	kyx=f3tensor(1,MM,1,NN,1,KK);
	kzx=f3tensor(1,MM,1,NN,1,KK);
	kxz=f3tensor(1,MM,1,NN,1,KK);
	kyz=f3tensor(1,MM,1,NN,1,KK);
	kzy=f3tensor(1,MM,1,NN,1,KK);
	spxx=matrix(1,MM,1,2*NN);
	spyy=matrix(1,MM,1,2*NN);
	spzz=matrix(1,MM,1,2*NN);
	spxy=matrix(1,MM,1,2*NN);
	spyx=matrix(1,MM,1,2*NN);
	spzx=matrix(1,MM,1,2*NN);
	spxz=matrix(1,MM,1,2*NN);
	spyz=matrix(1,MM,1,2*NN);
	spzy=matrix(1,MM,1,2*NN);
	
	for (i=1;i<=MM;i++){
		for (j=1;j<=NN;j++){
			for (k=1;k<=KK;k++){
				mm1[i][j][k]=0.0;
				mm2[i][j][k]=0.0;
				mm3[i][j][k]=0.0;
			}
		}
	}
	
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				mm1[i][j][k]=m1[i][j][k];
				mm2[i][j][k]=m2[i][j][k];
				mm3[i][j][k]=m3[i][j][k];
			}
		}
	}
	getkernel(kxx,kyy,kzz,kxy,kyz,kzx,spxx,spyy,spzz,spxy,spyz,spzx);
	for (i=1;i<=MM;i++){
		for (j=1;j<=NN;j++){
			for (k=1;k<=KK;k++){
	         kyx[i][j][k]=kxy[i][j][k];
	         kzy[i][j][k]=kyz[i][j][k];
	         kxz[i][j][k]=kzx[i][j][k];
	  		}
	  	}
	}
	for (i=1;i<=MM;i++){
		for (j=1;j<=2*NN;j++){
	         spyx[i][j]=spxy[i][j];
	         spzy[i][j]=spyz[i][j];
	         spxz[i][j]=spzx[i][j];
	  	}
	}	
	conv3(kxx,spxx,mm1,MM,NN,KK);
	conv3(kxy,spxy,mm2,MM,NN,KK);	
	conv3(kxz,spxz,mm3,MM,NN,KK);
	
	for (i=1;i<=MM;i++){
		for (j=1;j<=NN;j++){
			for (k=1;k<=KK;k++){
				mm1[i][j][k]=0.0;
				mm2[i][j][k]=0.0;
				mm3[i][j][k]=0.0;
			}
		}
	}
	
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				mm1[i][j][k]=m1[i][j][k];
				mm2[i][j][k]=m2[i][j][k];
				mm3[i][j][k]=m3[i][j][k];
			}
		}
	}
	
	conv3(kyx,spyx,mm1,MM,NN,KK);
	conv3(kyy,spyy,mm2,MM,NN,KK);
	conv3(kyz,spyz,mm3,MM,NN,KK);
	
   for (i=1;i<=MM;i++){
		for (j=1;j<=NN;j++){
			for (k=1;k<=KK;k++){
				mm1[i][j][k]=0.0;
				mm2[i][j][k]=0.0;
				mm3[i][j][k]=0.0;
			}
		}
	}
		
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				mm1[i][j][k]=m1[i][j][k];
				mm2[i][j][k]=m2[i][j][k];
				mm3[i][j][k]=m3[i][j][k];
			}
		}
	}
	
	conv3(kzx,spzx,mm1,MM,NN,KK);
	conv3(kzy,spzy,mm2,MM,NN,KK);
	conv3(kzz,spzz,mm3,MM,NN,KK);
	
	for(i=1;i<=M;i++){
		for(j=1;j<=N;j++){
			for(k=1;k<=K;k++){
				hd1[i][j][k]=mspi*(kxx[M+i][N+j][K+k]+
									  kxy[M+i][N+j][K+k]+
									  kxz[M+i][N+j][K+k]);
				hd2[i][j][k]=mspi*(kyx[M+i][N+j][K+k]+
									  kyy[M+i][N+j][K+k]+
									  kyz[M+i][N+j][K+k]);
				hd3[i][j][k]=mspi*(kzx[M+i][N+j][K+k]+
									  kzy[M+i][N+j][K+k]+
									  kzz[M+i][N+j][K+k]);
		   }
		}
	}
	
	free_f3tensor(mm1,1,MM,1,NN,1,KK);
	free_f3tensor(mm2,1,MM,1,NN,1,KK);
	free_f3tensor(mm3,1,MM,1,NN,1,KK);
	free_f3tensor(kxx,1,MM,1,NN,1,KK);
	free_f3tensor(kyy,1,MM,1,NN,1,KK);
	free_f3tensor(kzz,1,MM,1,NN,1,KK);
	free_f3tensor(kxy,1,MM,1,NN,1,KK);
	free_f3tensor(kyx,1,MM,1,NN,1,KK);
	free_f3tensor(kzx,1,MM,1,NN,1,KK);
	free_f3tensor(kxz,1,MM,1,NN,1,KK);
	free_f3tensor(kyz,1,MM,1,NN,1,KK);
	free_f3tensor(kzy,1,MM,1,NN,1,KK);
	free_matrix(spxx,1,MM,1,2*NN);
	free_matrix(spyy,1,MM,1,2*NN);
	free_matrix(spzz,1,MM,1,2*NN);
	free_matrix(spxy,1,MM,1,2*NN);
	free_matrix(spyx,1,MM,1,2*NN);
	free_matrix(spzx,1,MM,1,2*NN);
	free_matrix(spxz,1,MM,1,2*NN);
	free_matrix(spyz,1,MM,1,2*NN);
	free_matrix(spzy,1,MM,1,2*NN);
	
}	

#endif
	
	
