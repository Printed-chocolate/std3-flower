#ifndef EFFECTIVE_C
#define EFFECTIVE_C

#include "nrutil.c"

extern double gkms,gu,t,dx,dy,dz,mu_zero,Ms;
extern int M,N,K;

void effective(double ***m1,double ***m2,double ***m3,
	double ***h1,double ***h2,double ***h3,
	double ***he1,double ***he2,double ***he3,
	double ***f1,double ***f2,double ***f3,
	double coef)
{
	int i,j,k;
	double mx,my,mz,m12,m13,m21,m23,m31,m32;
	double alpha_m,alphax,alphay,alphaz;
	
	/***Here we can insert the externel field He which may
	depends on time t***/
	
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				mx=m1[i][j][k];
				my=m2[i][j][k];
				mz=m3[i][j][k];
                alphax=1;
				alphay=0;
				alphaz=0;
				alpha_m=-2.0*(mx*alphax+my*alphay+mz*alphaz);
				if(j<N)m12=m1[i][j+1][k];
				else m12=m1[i][N][k];
				if(k<K)m13=m1[i][j][k+1];
				else m13=m1[i][j][K];
				if(i<M)m21=m2[i+1][j][k];
				else m21=m2[M][j][k];
				if(k<K)m23=m2[i][j][k+1];
				else m23=m2[i][j][K];
				if(i<M)m31=m3[i+1][j][k];
				else m31=m3[M][j][k];
				if(j<N)m32=m3[i][j+1][k];
				else m32=m3[i][N][k];
				f1[i][j][k]=coef*(-gkms*(alpha_m*alphax) 
				/* mx*(my*my+mz*mz+2.0/4.6*my*my*mz*mz) */
					+gu*(h1[i][j][k]+2.0*he1[i][j][k])-2.0*0.005*(dy*(m32-mz)-dz*(m23-my)));
				f2[i][j][k]=coef*(-gkms*(alpha_m*alphay) 
				/* my*(mx*mx+mz*mz+2.0/4.6*mx*mx*mz*mz) */
					+gu*(h2[i][j][k]+2.0*he2[i][j][k])-2.0*0.005*(dz*(m13-mx)-dx*(m21-mx)));
				f3[i][j][k]=coef*(-gkms*(alpha_m*alphaz) 
				/* mz*(mx*mx+my*my+2.0/4.6*mx*mx*my*my) */
					+gu*(h3[i][j][k]+2.0*he3[i][j][k])-2.0*0.005*(dx*(m21-my)-dy*(m12-mx)));				
				
		   }
		}
	}
}

#endif		
