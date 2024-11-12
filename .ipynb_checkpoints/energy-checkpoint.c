#ifndef ENERGY_C
#define ENERGY_C

#include <stdlib.h>
#include <math.h>
#include "nrutil.c"

#define PI 3.141592653589793

extern double dx,dy,dz,Ms,A,Ku,mu_zero;
extern int M,N,K;

void energy(double *en, double *stray, double *exch,
	double ***m1,double ***m2,double ***m3,
	double ***hd1,double ***hd2,double ***hd3,
	double ***he1,double ***he2,double ***he3)
{
	int i,j,k;
        double alphax,alphay,alphaz;
	double ani,ani1,ani2,external,dv;
	double m10,m11,m20,m21,m30,m31,mx,my,mz;
	double ***ddm1,***ddm2,***ddm3;
	ani=0.0;ani1=0.0;ani2=0.0;*exch=0.0;*stray=0.0;external=0.0;
	dv=dx*dy*dz;
        alphax=1;
	alphay=0;
	alphaz=0;


	ddm1=f3tensor(1,M,1,N,1,K);
	ddm2=f3tensor(1,M,1,N,1,K);
	ddm3=f3tensor(1,M,1,N,1,K);
	
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				if(i>1)m10=m1[i-1][j][k];
				else m10=m1[1][j][k];
				if(i<M)m11=m1[i+1][j][k];
				else m11=m1[M][j][k];
				if(j>1)m20=m1[i][j-1][k];
				else m20=m1[i][1][k];
				if(j<N)m21=m1[i][j+1][k];
				else m21=m1[i][N][k];
				if(k>1)m30=m1[i][j][k-1];
				else m30=m1[i][j][1];
				if(k<K)m31=m1[i][j][k+1];
				else m31=m1[i][j][K];	
				ddm1[i][j][k]=pow((m11-m10)/2.0/dx,2)+
					pow((m21-m20)/2.0/dy,2)+
					pow((m31-m30)/2.0/dz,2);			
			}
		}
	}
   for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				if(i>1)m10=m2[i-1][j][k];
				else m10=m2[1][j][k];
				if(i<M)m11=m2[i+1][j][k];
				else m11=m2[M][j][k];
				if(j>1)m20=m2[i][j-1][k];
				else m20=m2[i][1][k];
				if(j<N)m21=m2[i][j+1][k];
				else m21=m2[i][N][k];
				if(k>1)m30=m2[i][j][k-1];
				else m30=m2[i][j][1];
				if(k<K)m31=m2[i][j][k+1];
				else m31=m2[i][j][K];	
				ddm2[i][j][k]=pow((m11-m10)/2.0/dx,2)+
					pow((m21-m20)/2.0/dy,2)+
					pow((m31-m30)/2.0/dz,2);	
				
			}
		}
	}
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				if(i>1)m10=m3[i-1][j][k];
				else m10=m3[1][j][k];
				if(i<M)m11=m3[i+1][j][k];
				else m11=m3[M][j][k];
				if(j>1)m20=m3[i][j-1][k];
				else m20=m3[i][1][k];
				if(j<N)m21=m3[i][j+1][k];
				else m21=m3[i][N][k];
				if(k>1)m30=m3[i][j][k-1];
				else m30=m3[i][j][1];
				if(k<K)m31=m3[i][j][k+1];
				else m31=m3[i][j][K];	
				ddm3[i][j][k]=pow((m11-m10)/2.0/dx,2)+
					pow((m21-m20)/2.0/dy,2)+
					pow((m31-m30)/2.0/dz,2);	
				
			}
		}
	}
	
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				mx=m1[i][j][k];my=m2[i][j][k];
				mz=m3[i][j][k];
			/*	ani1=dv*(mx*mx*my*my+my*my*mz*mz
					+mz*mz*mx*mx)
					+ani1; 
				ani2=dv*mx*mx*my*my*mz*mz
					+ani2; */    
			//	ani1=mx*alphax+my*alphay+mz*alphaz;
			//	ani2=ani1*ani1;
				ani=ani+dv*(1-mx*mx);
				*stray=*stray+
					dv*(hd1[i][j][k]*m1[i][j][k]+
					hd2[i][j][k]*m2[i][j][k]+
					hd3[i][j][k]*m3[i][j][k]);
				external=external+
					dv*(he1[i][j][k]*m1[i][j][k]+
					he2[i][j][k]*m2[i][j][k]+
					he3[i][j][k]*m3[i][j][k]);
				*exch=*exch+dv*(ddm1[i][j][k]+
					ddm2[i][j][k]+ddm3[i][j][k]);
			}
		}
	}
	*exch=A*(*exch);
	*stray=-0.5*mu_zero*Ms*(*stray);
//	*en=0.5*Ku*(ani1+1/2.3*ani2)+*exch+*stray
//		-mu_zero*Ms*external;
       	*en=Ku*(ani)+*exch+*stray
		-mu_zero*Ms*external;
	
	free_f3tensor(ddm1,1,M,1,N,1,K);
	free_f3tensor(ddm2,1,M,1,N,1,K);
	free_f3tensor(ddm3,1,M,1,N,1,K);
	
}	

#endif
	
	
