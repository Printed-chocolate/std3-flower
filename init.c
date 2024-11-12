#ifndef INIT_C
#define INIT_C

#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793

extern int M,N,K,new;
extern double dx,dy,dz,t,X,Y,Z,rad;
void init(double ***m1,double ***m2,double ***m3,
	double ***he1,double ***he2,double ***he3)
{
	FILE *fp;
	int i,j,k,idum1;
//	double G,H;
	double x,y,z,r,theta,norm,rand_max,a,aa;
    a=rad;
	aa=a*a;
	for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){
				for (k=1;k<=K;k++){
		            he1[i][j][k]=0.0;
		            he2[i][j][k]=0.0;
		            he3[i][j][k]=0.0;
		  		}
		  }
	}
	
	if(new==1){
		   t=0.0;
		   fp=fopen("x.dat","w");
		   if(fp==NULL){
				printf("FILE OPEN ERROR!\n");
				return;
		                }
		   for (i=1;i<=M;i++){
			fprintf(fp,"%e\n",(i-0.5)*dx);
		                      }
		   fclose(fp);
	
		   fp=fopen("y.dat","w");
	 	   if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		                }
		    for (j=1;j<=N;j++){
			fprintf(fp,"%e\n",(j-0.5)*dy);
		                      }
		   fclose(fp);
	
		   fp=fopen("z.dat","w");
		   if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		                }
		   for (k=1;k<=K;k++){
			fprintf(fp,"%e\n",(k-0.5)*dz);
		                     }
		   fclose(fp);
	
			for (j=1;j<=N;j++){
				for (k=1;k<=K;k++){
                    for (i=1;i<=M;i++){
						if(((i*dx-X/2)*(i*dx-X/2)+(j*dy-Y/2)*(j*dy-Y/2))<=aa){
							m1[i][j][k]=0.;
							m2[i][j][k]=0.;
							m3[i][j][k]=1.;
						} else {
							m1[i][j][k]=0.;
							m2[i][j][k]=0.;
							m3[i][j][k]=1.;
						}
						// if(((i*dx-140e-9)*(i*dx-140e-9)+(j*dy-64e-9)*(j*dy-64e-9))<=aa){
                        //     m1[i][j][k]=0.;
                        //     m2[i][j][k]=0.;
                        //     m3[i][j][k]=-1.;
                        // } else if(((i*dx-117e-9)*(i*dx-117e-9)+(j*dy-64e-9)*(j*dy-64e-9))<=aa){
                        //     m1[i][j][k]=0.;
                        //     m2[i][j][k]=0.;
                        //     m3[i][j][k]=-1.;
                        // } else if(((i*dx-94e-9)*(i*dx-94e-9)+(j*dy-64e-9)*(j*dy-64e-9))<=aa){
                        //     m1[i][j][k]=0.;
                        //     m2[i][j][k]=0.;
                        //     m3[i][j][k]=-1.;
                        // } else if(((i*dx-163e-9)*(i*dx-163e-9)+(j*dy-64e-9)*(j*dy-64e-9))<=aa){
                        //     m1[i][j][k]=0.;
                        //     m2[i][j][k]=0.;
                        //     m3[i][j][k]=-1.;
                        // } else {
                        //     m1[i][j][k]=0.;
                        //     m2[i][j][k]=0.;
                        //     m3[i][j][k]=1.;
                        // }	      
 						norm=sqrt(m1[i][j][k]*m1[i][j][k]+
				  			m2[i][j][k]*m2[i][j][k]+
				  			m3[i][j][k]*m3[i][j][k]);
						m1[i][j][k]=m1[i][j][k]/norm;
				  		m2[i][j][k]=m2[i][j][k]/norm;
						m3[i][j][k]=m3[i][j][k]/norm;
					}
				}
			}



		fp=fopen("m01.dat","w");
		   if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		            }
	  	   fprintf(fp,"%%t=%e\n",t);
		   for (k=1;k<=K;k++){
		       for (i=1;i<=M;i++){
			   for (j=1;j<=N;j++){
					fprintf(fp,"%e ",m1[i][j][k]);
				             }
				        fprintf(fp,"\n");
			                 }
		                      }
		fclose(fp);
	
		fp=fopen("m02.dat","w");
		   if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		                }
		fprintf(fp,"%%t=%e\n",t);
		for(k=1;k<=K;k++){	
		   for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){
					fprintf(fp,"%e ",m2[i][j][k]);
				           }
				fprintf(fp,"\n");
			              }
	 	                 }
		fclose(fp);
	
		fp=fopen("m03.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		             }
		fprintf(fp,"%%t=%e\n",t);
		for(k=1;k<=K;k++){
		    for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){
					fprintf(fp,"%e ",m3[i][j][k]);
				           }
				fprintf(fp,"\n");
			               }
		                  }
		fclose(fp);
		printf("The computation has been started from initial values!\n");
	   }		
          else{
   	        fp=fopen("m1.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		fscanf(fp,"%*c%*c%*c%le ",&t);
		for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){
				for (k=1;k<=K;k++){
					fscanf(fp,"%le ",&m1[i][j][k]);
				}
			}
		}
		fclose(fp);
	   fp=fopen("m2.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		fscanf(fp,"%*c%*c%*c%le ",&t);
		for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){
				for (k=1;k<=K;k++){	
					fscanf(fp,"%le ",&m2[i][j][k]);
				}
			}
		}
		fclose(fp);
		fp=fopen("m3.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		fscanf(fp,"%*c%*c%*c%le ",&t);
		for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){
				for (k=1;k<=K;k++){	
					fscanf(fp,"%le ",&m3[i][j][k]);
				}
			}
		}
		fclose(fp);
      printf("The computation has been started from last values!\n");
	}
}
#undef PI
#endif
