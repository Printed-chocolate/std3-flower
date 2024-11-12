#ifndef SOLVER_C
#define SOLVER_C
#include "nrutil.c"
#include "x_solver.c"
#include "xp_solver.c"
#include "y_solver.c"
#include "yp_solver.c"
#include "z_solver.c"
#include "xsolver.c"
#include "ysolver.c"
#include "zsolver.c"

/*Solve A*Lap(u)=F based on ADI; Results stored in u1.
 M,N,K are points number along X,Y,Z.
 Adtx=A*dt/3.0/dx/dx,etc.*/

//***solver1-Neumann boundary for XYZ (cosine transfer)

void solver1(double ***u, double ***u1,double ***F,
	unsigned long M, unsigned long N,unsigned long K,
	double Adtx,double Adty,double Adtz,double dt)
{
	double ***h;
	h=f3tensor(1,M,1,N,1,K);
	x_solver(u,u1,F,M,N,K,Adtx,Adty,Adtz,dt);
	y_solver(u1,h,F,M,N,K,Adtx,Adty,Adtz,dt);
	z_solver(h,u1,F,M,N,K,Adtx,Adty,Adtz,dt);
	free_f3tensor(h,1,M,1,N,1,K);
}

//***solver2-Neumann boundary for XZ,Perodic for Y

void solver2(double ***u, double ***u1,double ***F,
	unsigned long M, unsigned long N,unsigned long K,
	double Adtx,double Adty,double Adtz,double dt)
{
	double ***h;
	h=f3tensor(1,M,1,N,1,K);
	xp_solver(u,u1,F,M,N,K,Adtx,Adty,Adtz,dt);
	yp_solver(u1,h,F,M,N,K,Adtx,Adty,Adtz,dt);
	z_solver(h,u1,F,M,N,K,Adtx,Adty,Adtz,dt);
	free_f3tensor(h,1,M,1,N,1,K);
}

//***solver3-Neumann boundary for XYZ (Thomas method)

void solver3(double ***u, double ***u1,double ***F,
	unsigned long M, unsigned long N,unsigned long K,
	double Adtx,double Adty,double Adtz,double dt)
{
	double ***h;
	h=f3tensor(1,M,1,N,1,K);
	xsolver(u,u1,F,M,N,K,Adtx,Adty,Adtz,dt);
	ysolver(u1,h,F,M,N,K,Adtx,Adty,Adtz,dt);
	zsolver(h,u1,F,M,N,K,Adtx,Adty,Adtz,dt);
	free_f3tensor(h,1,M,1,N,1,K);
}


#endif
