#ifndef CONV3_C
#define CONV3_C
#include <stdlib.h>
#include "nrutil.c"
//#include "fourn.c"
#include "rlft3.c"
//This is the final version of 3-D convolution

void conv3(double ***data1,double **speq1,
	double ***data2,unsigned long n1,
	unsigned long n2,unsigned long n3)
{
   void rlft3(double ***data,double **speq, unsigned long nn1,
      unsigned long nn2, unsigned long nn3, int isign);
   int j;
   double fac,r,i,**speq2,*sp1,*sp2;

   //speq1=matrix(1,n1,1,2*n2);
   speq2=matrix(1,n1,1,2*n2);


   //rlft3(data1,speq1,n1,n2,n3,1);
   rlft3(data2,speq2,n1,n2,n3,1);
   fac=2.0/(n1*n2*n3);
   sp1=&data1[1][1][1];
   sp2=&data2[1][1][1];
   for (j=1;j<=n1*n2*n3/2;j++) {
   	r=sp1[0]*sp2[0]-sp1[1]*sp2[1];
   	i=sp1[0]*sp2[1]+sp1[1]*sp2[0];
   	sp1[0]=fac*r;
   	sp1[1]=fac*i;
   	sp1 +=2;
   	sp2+=2;
   }

   sp1=&speq1[1][1];
   sp2=&speq2[1][1];
   for (j=1;j<=n1*n2;j++) {
   	r=sp1[0]*sp2[0]-sp1[1]*sp2[1];
   	i=sp1[0]*sp2[1]+sp1[1]*sp2[0];
   	sp1[0]=fac*r;
   	sp1[1]=fac*i;
   	sp1+=2;
   	sp2+=2;
	}   	
 	rlft3(data1,speq1,n1,n2,n3,-1);
   	
   free_matrix(speq2,1,n1,1,2*n2);
  	//free_matrix(speq1,1,n1,1,2*n2);
  	  	
}

#endif
