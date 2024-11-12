#ifndef ENERGY_C
#define ENERGY_C

#include <stdlib.h>
#include <math.h>
#include "nrutil.c"

#define PI 3.141592653589793

extern double dx,dy,dz,Ms,A,Ku,mu_zero,gam,idmi;
extern int M,N,K;

void energy(double *en, double *stray, double *exch, double *ani, double *dmi,
	double ***m1,double ***m2,double ***m3,
	double ***hd1,double ***hd2,double ***hd3,
	double ***he1,double ***he2,double ***he3)
{
	int i,j,k;
    double alphax,alphay,alphaz;
	double ani1,ani2,external,dv;
	double m10,m11,m20,m21,m30,m31,mx,my,mz;
	double D, Dgu;
	double m1y,m1z,m3x,m3y,m2x,m2z;
	double ***ddm1,***ddm2,***ddm3;
	double ***hdmi1,***hdmi2,***hdmi3;
	*ani=0.0;ani1=0.0;ani2=0.0;*exch=0.0;*stray=0.0;external=0.0;*dmi=0.0;
	dv=dx*dy*dz;
    alphax=0;
	alphay=0;
	alphaz=1;
	ddm1=f3tensor(1,M,1,N,1,K);
	ddm2=f3tensor(1,M,1,N,1,K);
	ddm3=f3tensor(1,M,1,N,1,K);
	hdmi1=f3tensor(1,M,1,N,1,K);
	hdmi2=f3tensor(1,M,1,N,1,K);
	hdmi3=f3tensor(1,M,1,N,1,K);

	D = idmi;
	Dgu = -2.0*D*gam/Ms;
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				mx=m1[i][j][k];
				my=m2[i][j][k];
				mz=m3[i][j][k];
				double left_x = (j == 1) ? mx + D * dx * mz /(2*A) : m1[i][j-1][k];
				double right_x = (j == N) ? mx - D * dx * mz /(2*A) : m1[i][j+1][k];
				m1y = (right_x - left_x) / (2 * dy);
				double left1_x = (k == 1) ? mx - D * dx * my /(2*A) : m1[i][j][k-1];
				double right1_x = (k == M) ? mx + D * dx * my /(2*A) : m1[i][j][k+1];
				m1z = (right1_x - left1_x) / (2 * dz);
				double left_z = (i == 1) ? mz + D * dy * my /(2*A) : m3[i-1][j][k];
				double right_z = (i == M) ? mz - D * dy * my /(2*A) : m3[i+1][j][k];
				m3x = (right_z - left_z) / (2 * dx);
				double left1_z = (j == 1) ? mz - D * dx * mx /(2*A) : m3[i][j - 1][k];
				double right1_z = (j == N) ? mz + D * dx * mx /(2*A) : m3[i][j + 1][k];
				m3y = (right1_z - left1_z) / (2 * dy);
				double left_y = (i == 1) ? mz - D * dy * mz /(2*A) : m2[i-1][j][k];
				double right_y = (i == M) ? mz + D * dy * mz /(2*A) : m2[i+1][j][k];
				m2x = (right_y - left_y) / (2 * dx);
				double left1_y = (k == 1) ? my + D * dx * mx /(2*A) : m2[i][j][k-1];
				double right1_y = (k == K) ? my - D * dx * mx /(2*A) : m2[i][j][k+1];
				m2z = (right1_y - left1_y) / (2 * dz);
				hdmi1[i][j][k]=Dgu*(m2z-m3y);
				hdmi2[i][j][k]=Dgu*(m3x-m1z);
				hdmi3[i][j][k]=Dgu*(m1y-m2x);
		   }
		}
	}
	
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
				*ani=*ani+dv*(1-mz*mz);
				*dmi=*dmi+
                    dv*(hdmi1[i][j][k]*m1[i][j][k]+
					hdmi2[i][j][k]*m2[i][j][k]+
					hdmi3[i][j][k]*m3[i][j][k]);
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
	*ani=Ku*(*ani);
	*dmi=0.5*mu_zero*mu_zero*Ms*(*dmi);
	*stray=-0.5*mu_zero*Ms*(*stray);
    *en=*ani+*exch+*stray
		-mu_zero*Ms*external
		+*dmi;
	
	free_f3tensor(ddm1,1,M,1,N,1,K);
	free_f3tensor(ddm2,1,M,1,N,1,K);
	free_f3tensor(ddm3,1,M,1,N,1,K);
	free_f3tensor(hdmi1,1,M,1,N,1,K);
	free_f3tensor(hdmi2,1,M,1,N,1,K);
	free_f3tensor(hdmi3,1,M,1,N,1,K);
	
}	

#endif
	
	
