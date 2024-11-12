clear all;
close all;
load x.dat;
load y.dat;
load z.dat;
M=length(x);
N=length(y);
K=length(z);
red1=4;
red2=4;
red3=1;
xcomponent=length(x)/2;
ycomponent=length(y)/2;
zcomponent=length(z)/2;
MM1=load('mm1.dat');
MM2=load('mm2.dat');
MM3=load('mm3.dat');
%MM1=load('x1');
%MM2=load('y1');
%MM3=load('z1');

%load m01.dat;
%load m02.dat;
%load m03.dat;
%%[X,Y]=meshgrid(x,y);
n1=size(MM1,1);
n1=n1/K;
for k=1:K
    mm1(:,:,k)=MM1(((k-1)*n1+1):k*n1,:);
    mm2(:,:,k)=MM2(((k-1)*n1+1):k*n1,:);
    mm3(:,:,k)=MM3(((k-1)*n1+1):k*n1,:);
end
for i=1:M/red1 
   xr(i)=x(i*red1);
end
for j=1:N/red2
   yr(j)=y(j*red2);
end
for k=1:K/red3
   zr(k)=z(k*red3);
end
x=xr;
y=yr;
z=zr;
for i=1:M/red1
   for j=1:N/red2
      for k=1:K/red3
         mmr1(i,j,k)=mm1(i*red1,j*red2,k*red3);
         mmr2(i,j,k)=mm2(i*red1,j*red2,k*red3);
		mmr3(i,j,k)=mm3(i*red1,j*red2,k*red3);
      end
   end
end
mm1=mmr1;
mm2=mmr2;
mm3=mmr3;
xcomponent=xcomponent/red1
ycomponent=ycomponent/red2
zcomponent=zcomponent/red3
[x1,y1]=meshgrid(x,y);
[y2,z2]=meshgrid(y,z);
[x3,z3]=meshgrid(x,z);
m1=mm1(:,:,zcomponent);
m2=mm2(:,:,zcomponent);
m3=mm3(:,:,zcomponent);
%mx=zeros(N,K);
%my=zeros(N,K);
%mz=zeros(N,K);
%ma=zeros(M,K);
%mb=zeros(M,K);
%mc=zeros(M,K);
mx(:,:)=mm1(xcomponent,:,:);
my(:,:)=mm2(xcomponent,:,:);
mz(:,:)=mm3(xcomponent,:,:);
ma(:,:)=mm1(:,ycomponent,:);
mb(:,:)=mm2(:,ycomponent,:);
mc(:,:)=mm3(:,ycomponent,:);
figure;
%colormap('bone')
pcolor(x1,y1,m1');
shading interp;
%title('Top surface  M field , sample size 480x240x7.5 nm Arrows indicating (M1,M2),  Color indicating  M3');
%colorbar('vert');
axis('equal');
caxis([-1 1]);
hold    
quiver(x1,y1,m1',m2',0.5);
axis('equal');
xlabel('x')
ylabel('y')
figure;
%colormap('bone')
pcolor(y2,z2,mx');
shading interp;
colorbar('vert');
hold   
quiver(y2,z2,my',mz');
axis('equal');
%axis([0 4.8*10^(-7) -1*10^(-7) 1*10^(-7)])
caxis([-1 1]);
%hold
%figure;
%pcolor(x3,z3,mb');
%shading interp;
%colorbar('vert');
%hold
%quiver(x3,z3,ma',mc');
%axis('equal');
%figure;
%mesh(x1,y1,m3');
%%[x3,y3,z3]=meshgrid(x,y,z);
%%quiver3(x3,y3,z3,mm1,mm2,mm3);
%
%%%%%quiver(X,Y,m2,m1);
