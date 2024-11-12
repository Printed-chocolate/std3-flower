#ifndef OUTPUT_C
#define OUTPUT_C

#include <stdio.h>

int M,N,K;
double t;

void output(double ***m1,double ***m2,double ***m3,
		double en, double stray, double exch, double ani, double dmi)
{
	FILE *fp;
	int i,j,k;
	fp=fopen("view","a");
	if(fp==NULL){
		printf("FILE OPEN ERROR!\n");
		return;
	}
	fprintf(fp,"Time=%e ; Total Energy=%e ; Exchange Energy=%e; StrayField Energy=%e; Anisotropy Energy=%e; DMI Energy=%e\n",t,en,exch,stray,ani,dmi);  	
	fclose(fp);
		
	fp=fopen("m1.dat","w");
	if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
	}
	fprintf(fp,"%%t=%e\n",t);
	for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){
				for(k=1;k<=K;k++){
					fprintf(fp,"%.8e ",m1[i][j][k]);
			   }
				fprintf(fp,"\n");
			}
	}
	fclose(fp);
	
	fp=fopen("m2.dat","w");
	if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
	}
	fprintf(fp,"%%t=%e\n",t);
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for(k=1;k<=K;k++){
				fprintf(fp,"%.8e ",m2[i][j][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);
	
	fp=fopen("m3.dat","w");
	if(fp==NULL){
			printf("FILE OPEN ERROR!\n");
			return;
	}
	fprintf(fp,"%%t=%e\n",t);
	for (i=1;i<=M;i++){
		for (j=1;j<=N;j++){
			for(k=1;k<=K;k++){
				fprintf(fp,"%.8e ",m3[i][j][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);		
	
	fp=fopen("mm1.dat","w");
	if(fp==NULL){
		printf("FILE OPEN ERROR!\n");
		return;
	}
	for (k=1;k<=K;k++){
		for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){		
				fprintf(fp,"%e ",m1[i][j][k]);
			}
		   fprintf(fp,"\n");
		}
	}
	fclose(fp);
	fp=fopen("mm2.dat","w");
	if(fp==NULL){
		printf("FILE OPEN ERROR!\n");
		return;
	}
	for (k=1;k<=K;k++){
		for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){		
				fprintf(fp,"%e ",m2[i][j][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);

	fp=fopen("mm3.dat","w");
	if(fp==NULL){
		printf("FILE OPEN ERROR!\n");
		return;
	}
	for (k=1;k<=K;k++){
		for (i=1;i<=M;i++){
			for (j=1;j<=N;j++){		
				fprintf(fp,"%e ",m3[i][j][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fclose(fp);		
}
#endif	
