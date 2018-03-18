%Seth Reyes
clear all, clc

Lx= 1 ;
Ly= 1 ;
u0= 0 ;
uL= 1 ;
v0= 0 ;
vL= 1 ;
h=100; g=h;
dx=Lx/(h+1);
dy=Ly/(g+1);

%M=2;
%f = -2*M*sin(M*x)*cosh(M*y);
for i=2:h
    
    uMat(i)=1;
end