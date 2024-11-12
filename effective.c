#ifndef EFFECTIVE_C
#define EFFECTIVE_C

#include "nrutil.c"

extern double gkms,gu,t,dx,dy,dz,A,gam,Ms,mu_zero,dcd,sha,idmi;
extern int M,N,K;

void effective(double ***m1,double ***m2,double ***m3,
	double ***h1,double ***h2,double ***h3,
	double ***he1,double ***he2,double ***he3,
	double ***f1,double ***f2,double ***f3,
	double coef)
{
	int i,j,k;
	double mx,my,mz,m1y,m1z,m3x,m3y,m2x,m2z;
	double m1x,m2y;
	double alpha_m,alphax,alphay,alphaz,Dgu,qe,rp,sotcoefficient;
	
	/***Here we can insert the externel field He which may depends on time t***/
	alphax=0;
	alphay=0;
	alphaz=1;
	Dgu = -2*idmi*gam/Ms;
    rp = 1.055e-34;
    qe = 1.6e-19;
    sotcoefficient = gam*rp*dcd*sha/(2*qe*Ms*mu_zero*dz);
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for (k=1;k<=K;k++){
				mx=m1[i][j][k];
				my=m2[i][j][k];
				mz=m3[i][j][k];
				double left_x = (i == 1) ? mz - idmi * dx * mx /(2*A) : m3[i - 1][j][k];
				double right_x = (i == M) ? mz + idmi * dx * mx /(2*A) : m3[i + 1][j][k];
				m3x = (right_x - left_x) / (2 * dx);
				double left1_x = (i == 1) ? mx + idmi * dx * mz /(2*A) : m1[i - 1][j][k];
				double right1_x = (i == M) ? mx - idmi * dx * mz /(2*A) : m1[i + 1][j][k];
				m1x = (right1_x - left1_x) / (2 * dx);
				double left_y = (j == 1) ? mz - idmi * dy * my /(2*A) : m3[i][j - 1][k];
				double right_y = (j == N) ? mz + idmi * dy * my /(2*A) : m3[i][j + 1][k];
				m3y = (right_y - left_y) / (2 * dy);
				double left1_y = (j == 1) ? my + idmi * dx * mz /(2*A) : m2[i][j - 1][k];
				double right1_y = (j == N) ? my - idmi * dx * mz /(2*A) : m2[i][j + 1][k]; 
				m2y = (right1_y - left1_y) / (2 * dy);
				alpha_m=-2.0*(mx*alphax+my*alphay+mz*alphaz);
/*				f1[i][j][k]=coef*(-gkms*(alpha_m*alphax) 
					+gu*(h1[i][j][k]+2.0*he1[i][j][k])-gam*eta*jr*mz);
				f2[i][j][k]=coef*(-gkms*(alpha_m*alphay) 
					+gu*(h2[i][j][k]+2.0*he2[i][j][k]));
				f3[i][j][k]=coef*(-gkms*(alpha_m*alphaz) 
					+gu*(h3[i][j][k]+2.0*he3[i][j][k])+gam*eta*jr*mx+gam*zeta*e*mz);*/
				f1[i][j][k]=coef*(-gkms*(alpha_m*alphax) 
					+gu*(h1[i][j][k]+2.0*he1[i][j][k])+Dgu*m3x+sotcoefficient*mz);
				f2[i][j][k]=coef*(-gkms*(alpha_m*alphay) 
					+gu*(h2[i][j][k]+2.0*he2[i][j][k])+Dgu*m3y);
				f3[i][j][k]=coef*(-gkms*(alpha_m*alphaz) 
					+gu*(h3[i][j][k]+2.0*he3[i][j][k])-Dgu*m1x-Dgu*m2y-sotcoefficient*mx);
/*				double left_x = (j == 1) ? mx + D * dx * mz /(2*A) : m1[i][j-1][k];
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
				alpha_m=-2.0*(mx*alphax+my*alphay+mz*alphaz);
				f1[i][j][k]=coef*(-gkms*(alpha_m*alphax) 
					+gu*(h1[i][j][k]+2.0*he1[i][j][k])+Dgu*(m2z-m3y));
				f2[i][j][k]=coef*(-gkms*(alpha_m*alphay) 
					+gu*(h2[i][j][k]+2.0*he2[i][j][k])+Dgu*(m3x-m1z));
				f3[i][j][k]=coef*(-gkms*(alpha_m*alphaz) 
					+gu*(h3[i][j][k]+2.0*he3[i][j][k])+Dgu*(m1y-m2x)); */
		   }
		}
	}
} 

#endif

/* 				double left_x = (i == 1) ? mz - D * dx * mx /(2*A) : m3[i - 1][j][k];
				double right_x = (i == M) ? mz + D * dx * mx /(2*A) : m3[i + 1][j][k];
				m3x = (right_x - left_x) / (2 * dx);
				double left1_x = (i == 1) ? mx + D * dx * mz /(2*A) : m1[i - 1][j][k];
				double right1_x = (i == M) ? mx - D * dx * mz /(2*A) : m1[i + 1][j][k];
				m1x = (right1_x - left1_x) / (2 * dx);
				double left_y = (j == 1) ? mz - D * dy * my /(2*A) : m3[i][j - 1][k];
				double right_y = (j == N) ? mz + D * dy * my /(2*A) : m3[i][j + 1][k];
				m3y = (right_y - left_y) / (2 * dy);
				double left1_y = (j == 1) ? my + D * dx * mz /(2*A) : m2[i][j - 1][k];
				double right1_y = (j == N) ? my - D * dx * mz /(2*A) : m2[i][j + 1][k];
				m2y = (right1_y - left1_y) / (2 * dy);
				alpha_m=-2.0*(mx*alphax+my*alphay+mz*alphaz);
				f1[i][j][k]=coef*(-gkms*(alpha_m*alphax) 
					+gu*(h1[i][j][k]+2.0*he1[i][j][k])+Dgu*m3x);
				f2[i][j][k]=coef*(-gkms*(alpha_m*alphay) 
					+gu*(h2[i][j][k]+2.0*he2[i][j][k])+Dgu*m3y);
				f3[i][j][k]=coef*(-gkms*(alpha_m*alphaz) 
					+gu*(h3[i][j][k]+2.0*he3[i][j][k])-Dgu*m1x-Dgu*m2y); */