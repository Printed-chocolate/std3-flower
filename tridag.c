/* note #undef's at end of file */
#ifndef TRIDAG_C
#define TRIDAG_C
#define NRANSI
#include "nrutil.c"

void tridag(double a, double b, double c,
	double b1, double *r, double *u,
	unsigned long n)
{
	unsigned long j;
	double bet;
	double *gam;
	gam=vector(1,n);
	if (b1 == 0.0) nrerror("Error 1 in tridag");
	u[1]=r[1]/(bet=b1);
	for (j=2;j<=n-1;j++) {
	   gam[j]=c/bet;
	   bet=b-a*gam[j];
	   if (bet == 0.0)	nrerror("Error 2 in tridag");
		u[j]=(r[j]-a*u[j-1])/bet;
	}
	gam[n]=c/bet;
	bet=b1-a*gam[n];
	if (bet == 0.0)	nrerror("Error 2 in tridag");
	u[n]=(r[n]-a*u[n-1])/bet;
	for (j=(n-1);j>=1;j--)
		u[j] -= gam[j+1]*u[j+1];
}
	
#undef NRANSI
#endif
