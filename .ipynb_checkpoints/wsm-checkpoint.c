//This wsm.c is the final version and gsp.c is not!!!!!
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "nrutil.c"
#include "readdata.c"
#include "init.c"
#include "kernel.c"
#include "hd.c"
#include "effective.c"
#include "solver.c"
#include "output.c"
#include "energy.c"

#define PI 3.141592653589793


int M,N,K,MM,NN,KK;
int nsum,ker,new,boundary,npoints;
double X,Y,Z,dx,dy,dz,t,dt,nsteps,epsilon,RelTol;
double mu_zero,Ms,A,Ku,gam,alpha,gkms,gams,gu,mspi;




int main(void)
{
	double ***m1,***m2,***m3;
	double ***f1,***f2,***f3;
	double ***g1,***g2,***g3;
	double ***h1,***h2,***h3;
	double ***he1,***he2,***he3;
	
	void (*solver)();
	
	double Adtx,Adty,Adtz,Bdtx,Bdty,Bdtz;
	double norm,en,stray,exch,en0,ani;

    double unit,unitot,ustray,uniexch,uniani;
	double Kd;
	int i,j,k,count,reach;
	char str[80];
	FILE *fp,*fp1;
	
	//***read data from file data.c
	
	readdata();
	
	//***Initialize some parameters
	
   MM=2*M;
   NN=2*N;
   KK=2*K;
   dx=X/(M+0.0);
   dy=Y/(N+0.0);
   dz=Z/(K+0.0);
	
   mu_zero=4.0*PI*1.0e-7;
   gkms=gam*Ku/Ms;
   gams=gam*A/Ms;
   gu=0.5*gam*mu_zero;
   mspi=Ms/4.0/PI;
   Kd=mu_zero*Ms*Ms/2;

   Adtx=gams*dt/dx/dx/3.0;
   Adty=gams*dt/dy/dy/3.0;
   Adtz=gams*dt/dz/dz/3.0;
   Bdtx=Adtx*alpha;
   Bdty=Adty*alpha;
   Bdtz=Adtz*alpha;

   //***Magnetization M=(m1,m2,m3)

        m1=f3tensor(1,M,1,N,1,K);
	m2=f3tensor(1,M,1,N,1,K);
	m3=f3tensor(1,M,1,N,1,K);
	
   //***Stray field Hs=(h1,h2,h3)
	
        h1=f3tensor(1,M,1,N,1,K);
	h2=f3tensor(1,M,1,N,1,K);
	h3=f3tensor(1,M,1,N,1,K);
	
	//***External field He=(he1,he2,he3)
	
	he1=f3tensor(1,M,1,N,1,K);
	he2=f3tensor(1,M,1,N,1,K);
	he3=f3tensor(1,M,1,N,1,K);

        f1=f3tensor(1,M,1,N,1,K);
	f2=f3tensor(1,M,1,N,1,K);
	f3=f3tensor(1,M,1,N,1,K);	
	
	g1=f3tensor(1,M,1,N,1,K);
	g2=f3tensor(1,M,1,N,1,K);
	g3=f3tensor(1,M,1,N,1,K);
		
	//***Generating or loading kernels forstray field
	
	kernel(nsum);
	
	//***Initializing magnetiztion and external field
	
	init(m1,m2,m3,he1,he2,he3);
	
	//***Choose Neumann boundary (cosine transformation)
	
	if (boundary==1)
		solver=solver1;
	
	//***X&Z-Neumann boundary; Y-Periodic	
	
	else if (boundary==2)
		solver=solver2;
	
	//***Choose Neumann boundary (Thomas method)	
	
	else if (boundary==3)
		solver=solver3;
	
	//***The main loop of time advances by using G-S-P method
	
	for (count=1;count<=nsteps;count++){
		t=t+1.0*dt;
		
		//**Compute stray field
		
		hd(m1,m2,m3,h1,h2,h3);
		
		//**Compute effective field (with no exchange energy)
		
		effective(m1,m2,m3,h1,h2,h3,
			he1,he2,he3,f1,f2,f3,1.0);
		
		/**Solve heat equation (with source) subject to chosen
		 boundary condition**/
		
		(*solver)(m2,g2,f2,M,N,K,Adtx,Adty,Adtz,dt);
		(*solver)(m3,g3,f3,M,N,K,Adtx,Adty,Adtz,dt);
		
		
		//**Update m1 by G-S iteration
		
		for(i=1;i<=M;i++){
			for(j=1;j<=N;j++){
				for(k=1;k<=K;k++){
					m1[i][j][k]=m1[i][j][k]
						+(g2[i][j][k]*m3[i][j][k]
						-g3[i][j][k]*m2[i][j][k]);
			   }
			}
		}
		
		//**Update effective field
		
		effective(m1,m2,m3,h1,h2,h3,
			he1,he2,he3,f1,f2,f3,1.0);
		
		
		(*solver)(m1,g1,f1,M,N,K,Adtx,Adty,Adtz,dt);
		
		//**Update m2 by G-S iteration
		
		for(i=1;i<=M;i++){
			for(j=1;j<=N;j++){
				for(k=1;k<=K;k++){
					m2[i][j][k]=m2[i][j][k]
						+(g3[i][j][k]*m1[i][j][k]
						-g1[i][j][k]*m3[i][j][k]);
				}
			}
		}
		
      effective(m1,m2,m3,h1,h2,h3,
      	he1,he2,he3,f1,f2,f3,1.0);
				
		(*solver)(m2,g2,f2,M,N,K,Adtx,Adty,Adtz,dt);
		
		//**Update m3 by G-S iteration
		
		for(i=1;i<=M;i++){
			for(j=1;j<=N;j++){
				for(k=1;k<=K;k++){		
					m3[i][j][k]=m3[i][j][k]
						+(g1[i][j][k]*m2[i][j][k]
						-g2[i][j][k]*m1[i][j][k]);
				}
			}
		}
		
		//**Solve heat equation without constrain
		
		effective(m1,m2,m3,h1,h2,h3,
			he1,he2,he3,f1,f2,f3,alpha);
		
		(*solver)(m1,g1,f1,M,N,K,Bdtx,Bdty,Bdtz,dt);
		(*solver)(m2,g2,f2,M,N,K,Bdtx,Bdty,Bdtz,dt);
		(*solver)(m3,g3,f3,M,N,K,Bdtx,Bdty,Bdtz,dt);
		
		//**Project to unit sphere
		
		for(i=1;i<=M;i++){
			for(j=1;j<=N;j++){
				for(k=1;k<=K;k++){
					norm=sqrt(g1[i][j][k]*g1[i][j][k]
						+g2[i][j][k]*g2[i][j][k]
						+g3[i][j][k]*g3[i][j][k]);
					m1[i][j][k]=g1[i][j][k]/norm;
					m2[i][j][k]=g2[i][j][k]/norm;
					m3[i][j][k]=g3[i][j][k]/norm;
			
			   }
			}
		}
	
	   //**Output m1,m2,m3 and total energy

	        
	
	
	   if(fmod(count,nsteps/(npoints+0.0))==0.0){
		   en0=en;
	   	   energy(&en,&stray,&exch,m1,m2,m3,h1,h2,h3,he1,he2,he3);
                   output(m1,m2,m3,en,stray,exch,ani);
		   if(fabs((en-en0)/en0)<=RelTol)
		       {
	   	         break;
                         output(m1,m2,m3,en,stray,exch,ani);
		       }
                  
	   }
	}
	
	//***Mark of termination of program
	        
                unit=Kd*X*Y*Z;
		unitot=en/unit;
		ustray=stray/unit;
		uniexch=exch/unit;
		uniani=unitot-ustray-uniexch;
	fp=fopen("view","a");
	fprintf(fp,"========END========\n");
	fprintf(fp," Unit Energy=%e ;\n StrayField Energy=%e ;\n Exchange Energy=%e ;\n Anisotropy Energy=%e\n",unitot,ustray,uniexch,uniani);
	fclose(fp);
	
	//***Free kernels from memory
	
	freekernel();
	
	//***Free memories
	
	free_f3tensor(m1,1,M,1,N,1,K);
	free_f3tensor(m2,1,M,1,N,1,K);
	free_f3tensor(m3,1,M,1,N,1,K);
	free_f3tensor(f1,1,M,1,N,1,K);
	free_f3tensor(f2,1,M,1,N,1,K);
	free_f3tensor(f3,1,M,1,N,1,K);
	free_f3tensor(g1,1,M,1,N,1,K);
	free_f3tensor(g2,1,M,1,N,1,K);
	free_f3tensor(g3,1,M,1,N,1,K);
	free_f3tensor(h1,1,M,1,N,1,K);
	free_f3tensor(h2,1,M,1,N,1,K);
	free_f3tensor(h3,1,M,1,N,1,K);
	free_f3tensor(he1,1,M,1,N,1,K);
	free_f3tensor(he2,1,M,1,N,1,K);
	free_f3tensor(he3,1,M,1,N,1,K);
	
	
}
#undef PI
	
	
