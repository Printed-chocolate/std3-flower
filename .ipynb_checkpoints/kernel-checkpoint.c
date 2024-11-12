#ifndef KERNEL_C
#define KERNEL_C	

#include <math.h>
#include "nrutil.c"
#include "rlft3.c"

static double ***Kxx, ***Kyy,***Kzz,***Kxy,***Kyz,***Kzx;
static double **Spxx,**Spyy,**Spzz,**Spxy,**Spyz,**Spzx;
extern double Y,dx,dy,dz;
extern int M,N,K,MM,NN,KK;
extern int ker;


double fx(double x, double y, double z)
{
	return atan(y*z/x/sqrt(x*x+y*y+z*z));	
}

double fy(double x, double y, double z)
{
	return atan(x*z/y/sqrt(x*x+y*y+z*z));	
}

double fz(double x, double y, double z)
{
   return atan(x*y/z/sqrt(x*x+y*y+z*z));	
}

double gz(double x,double y,double z)
{
	return -log(z+sqrt(x*x+y*y+z*z));
}

double gx(double x,double y,double z)
{
	return -log(x+sqrt(x*x+y*y+z*z));
}

double gy(double x,double y,double z)
{
	return -log(y+sqrt(x*x+y*y+z*z));
}

double Fxx(double x,double y1,double y2,
		double z1,double z2)
{
	return  fx(x,y2,z2)-fx(x,y2,z1)
			-fx(x,y1,z2)+fx(x,y1,z1);
}

double Fyy(double y,double x1,double x2,
		double z1,double z2)
{
	return  fy(x2,y,z2)-fy(x2,y,z1)
		-fy(x1,y,z2)+fy(x1,y,z1);
}

double Fzz(double z,double x1,double x2,
		double y1,double y2)
{
	return  fz(x2,y2,z)-fz(x1,y2,z)
		-fz(x2,y1,z)+fz(x1,y1,z);
}

double Fyz(double z,double x1,double x2,
		double y1,double y2)
{
	return  gx(x2,y2,z)-gx(x2,y1,z)
		-gx(x1,y2,z)+gx(x1,y1,z);
}

double Fzx(double x,double y1,double y2,
		double z1,double z2)
{
	return  gy(x,y2,z2)-gy(x,y2,z1)
		-gy(x,y1,z2)+gy(x,y1,z1);
}

double Fxy(double y,double x1,double x2,
		double z1,double z2)
{
	return  gz(x2,y,z2)-gz(x2,y,z1)
		-gz(x1,y,z2)+gz(x1,y,z1);
}



void kernel(int nsum)
{
	int p,q,r,n;
	double x1,x2,y1,y2,z1,z2;
	FILE *fp;
	Kxx=f3tensor(1,MM,1,NN,1,KK);
	Kyy=f3tensor(1,MM,1,NN,1,KK);
	Kzz=f3tensor(1,MM,1,NN,1,KK);
	Kxy=f3tensor(1,MM,1,NN,1,KK);
	Kyz=f3tensor(1,MM,1,NN,1,KK);
	Kzx=f3tensor(1,MM,1,NN,1,KK);
	Spxx=matrix(1,MM,1,2*NN);
	Spyy=matrix(1,MM,1,2*NN);
	Spzz=matrix(1,MM,1,2*NN);
	Spxy=matrix(1,MM,1,2*NN);
	Spyz=matrix(1,MM,1,2*NN);
	Spzx=matrix(1,MM,1,2*NN);
	
	if (ker==1) {
		printf("The kernel is being calculated for the first time!\n");
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for (r=1;r<=KK;r++){
					Kxx[p][q][r]=0.0;
					Kyy[p][q][r]=0.0;
					Kzz[p][q][r]=0.0;
					Kxy[p][q][r]=0.0;
					Kyz[p][q][r]=0.0;
					Kzx[p][q][r]=0.0;
			   }
			}
		}
		for (n=-nsum;n<=nsum;n++){
			for(p=1-M;p<=M-1;p++){
				for (q=1-N;q<=N-1;q++){
					for(r=1-K;r<=K-1;r++){
						x1=dx*(p+0.5);x2=x1-dx;
						y1=dy*(q+0.5)-n*Y;y2=y1-dy;
						z1=dz*(r+0.5);z2=z1-dz;
						Kxx[p+M+1][q+N+1][r+K+1]=
							Fxx(x2,y1,y2,z1,z2)
							-Fxx(x1,y1,y2,z1,z2)
							+Kxx[p+M+1][q+N+1][r+K+1];
						Kyy[p+M+1][q+N+1][r+K+1]=
							Fyy(y2,x1,x2,z1,z2)
							-Fyy(y1,x1,x2,z1,z2)
							+Kyy[p+M+1][q+N+1][r+K+1];
						Kzz[p+M+1][q+N+1][r+K+1]=
							Fzz(z2,x1,x2,y1,y2)
							-Fzz(z1,x1,x2,y1,y2)
							+Kzz[p+M+1][q+N+1][r+K+1];
						Kxy[p+M+1][q+N+1][r+K+1]=
							Fxy(y2,x1,x2,z1,z2)
							-Fxy(y1,x1,x2,z1,z2)
							+Kxy[p+M+1][q+N+1][r+K+1];
						Kyz[p+M+1][q+N+1][r+K+1]=
							Fyz(z2,x1,x2,y1,y2)
							-Fyz(z1,x1,x2,y1,y2)
							+Kyz[p+M+1][q+N+1][r+K+1];
						Kzx[p+M+1][q+N+1][r+K+1]=
							Fzx(x2,y1,y2,z1,z2)
							-Fzx(x1,y1,y2,z1,z2)
							+Kzx[p+M+1][q+N+1][r+K+1];
					}
				}
			}
		}	
   	
   	rlft3(Kxx,Spxx,MM,NN,KK,1);
   	rlft3(Kyy,Spyy,MM,NN,KK,1);
   	rlft3(Kzz,Spzz,MM,NN,KK,1);
   	rlft3(Kxy,Spxy,MM,NN,KK,1);
   	rlft3(Kyz,Spyz,MM,NN,KK,1);
   	rlft3(Kzx,Spzx,MM,NN,KK,1);
   				
/*		fp=fopen("Kxx.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				 for(r=1;r<=KK;r++){
					fprintf(fp,"%.8e ",Kxx[p][q][r]);
				 }
				 fprintf(fp,"\n");
			}
		}
		fclose(fp);
	   fp=fopen("Kyy.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				 for(r=1;r<=KK;r++){
					fprintf(fp,"%.8e ",Kyy[p][q][r]);
				 }
				 fprintf(fp,"\n");
			}
		}
		fclose(fp);
		fp=fopen("Kzz.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				 for(r=1;r<=KK;r++){
					fprintf(fp,"%.8e ",Kzz[p][q][r]);
				 }
				 fprintf(fp,"\n");
			}
		}
		fclose(fp);
		fp=fopen("Kxy.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				 for(r=1;r<=KK;r++){
					fprintf(fp,"%.8e ",Kxy[p][q][r]);
				 }
				 fprintf(fp,"\n");
			}
		}
		fclose(fp);
		fp=fopen("Kyz.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for(r=1;r<=KK;r++){
					fprintf(fp,"%.8e ",Kyz[p][q][r]);
				 }
				 fprintf(fp,"\n");
			}
		}
		fclose(fp);
		fp=fopen("Kzx.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for(r=1;r<=KK;r++){
					fprintf(fp,"%.8e ",Kzx[p][q][r]);
				 }
				 fprintf(fp,"\n");
			}
		}
		fclose(fp);
		fp=fopen("spxx.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fprintf(fp,"%.8e ",Spxx[p][q]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		fp=fopen("spyy.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fprintf(fp,"%.8e ",Spyy[p][q]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		fp=fopen("spzz.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fprintf(fp,"%.8e ",Spzz[p][q]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		fp=fopen("spxy.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fprintf(fp,"%.8e ",Spxy[p][q]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		fp=fopen("spyz.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fprintf(fp,"%.8e ",Spyz[p][q]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		fp=fopen("spzx.dat","w");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fprintf(fp,"%.8e ",Spzx[p][q]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);*/
	}
	else{
		printf("The kernel is loaded from the files!\n");
/*		fp=fopen("Kxx.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for(r=1;r<=KK;r++){	
					fscanf(fp,"%le ",&Kxx[p][q][r]);
				}
			}
		}
		fclose(fp);
	   fp=fopen("Kyy.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for(r=1;r<=KK;r++){	
					fscanf(fp,"%le ",&Kyy[p][q][r]);
				}
			}
		}
		fclose(fp);
		fp=fopen("Kzz.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for(r=1;r<=KK;r++){	
					fscanf(fp,"%le ",&Kzz[p][q][r]);
				}
			}
		}
		fclose(fp);
		fp=fopen("Kxy.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for(r=1;r<=KK;r++){	
					fscanf(fp,"%le ",&Kxy[p][q][r]);
				}
			}
		}
		fclose(fp);
		fp=fopen("Kyz.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for(r=1;r<=KK;r++){	
					fscanf(fp,"%le ",&Kyz[p][q][r]);
				}
			}
		}
		fclose(fp);
		fp=fopen("Kzx.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=NN;q++){
				for(r=1;r<=KK;r++){	
					fscanf(fp,"%le ",&Kzx[p][q][r]);
				}
			}
		}
		fclose(fp);
		fp=fopen("spxx.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fscanf(fp,"%le ",&Spxx[p][q]);
			}
		}
		fclose(fp);
		fp=fopen("spyy.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fscanf(fp,"%le ",&Spyy[p][q]);
			}
		}
		fclose(fp);
		fp=fopen("spzz.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fscanf(fp,"%le ",&Spzz[p][q]);
			}
		}
		fclose(fp);
		fp=fopen("spxy.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fscanf(fp,"%le ",&Spxy[p][q]);
			}
		}
		fclose(fp);
		fp=fopen("spyz.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fscanf(fp,"%le ",&Spyz[p][q]);
			}
		}
		fclose(fp);
		fp=fopen("spzx.dat","r");
		if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
		}
		for(p=1;p<=MM;p++){
			for(q=1;q<=2*NN;q++){
				fscanf(fp,"%le ",&Spzx[p][q]);
			}
		}
		fclose(fp);*/
	}
}

void getkernel(double ***kxx,double ***kyy,double ***kzz,
	double ***kxy,double ***kyz,double ***kzx,
	double **spxx,double **spyy,double **spzz,
	double **spxy,double **spyz,double **spzx)
{
	int i,j,k;
	for(i=1;i<=MM;i++){
		for(j=1;j<=NN;j++){
			for(k=1;k<=KK;k++){
				kxx[i][j][k]=Kxx[i][j][k];
				kyy[i][j][k]=Kyy[i][j][k];
				kzz[i][j][k]=Kzz[i][j][k];
				kxy[i][j][k]=Kxy[i][j][k];
				kyz[i][j][k]=Kyz[i][j][k];
				kzx[i][j][k]=Kzx[i][j][k];		
			}
		}
	}
	for(i=1;i<=MM;i++){
		for(j=1;j<=2*NN;j++){
			spxx[i][j]=Spxx[i][j];	
			spyy[i][j]=Spyy[i][j];
			spzz[i][j]=Spzz[i][j];
			spxy[i][j]=Spxy[i][j];
			spyz[i][j]=Spyz[i][j];
			spzx[i][j]=Spzx[i][j];
		}
	}
	
}

void freekernel(void)
{
	free_f3tensor(Kxx,1,MM,1,NN,1,KK);
	free_f3tensor(Kyy,1,MM,1,NN,1,KK);
	free_f3tensor(Kzz,1,MM,1,NN,1,KK);	
	free_f3tensor(Kxy,1,MM,1,NN,1,KK);
	free_f3tensor(Kyz,1,MM,1,NN,1,KK);
	free_f3tensor(Kzx,1,MM,1,NN,1,KK);	
	free_matrix(Spxx,1,MM,1,2*NN);
	free_matrix(Spyy,1,MM,1,2*NN);
	free_matrix(Spzz,1,MM,1,2*NN);
	free_matrix(Spxy,1,MM,1,2*NN);
	free_matrix(Spyz,1,MM,1,2*NN);
	free_matrix(Spzx,1,MM,1,2*NN);
}

#endif
	
