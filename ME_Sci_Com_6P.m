%Seth Reyes
clear all, clc,close all

Lx= pi ;
Ly= pi ;
u0= 0 ;
uL= 0 ;
v0= 0 ;
vL= 0 ;
h=25; g=h;
dx=Lx/(h);
dy=Ly/(g);
M=1;
Adiag=eye(h-1)*(-4);
Adiag=Adiag+diag(ones(h-2,1),1);
Adiag=Adiag+diag(ones(h-2,1),-1);
Aother=eye(h-1);
A=zeros((h-1)*(h-1),(h-1)*(h-1));
for i=1:h-1
    A((i-1)*(h-1)+1:(i-1)*(h-1)+(h-1),(i-1)*(h-1)+1:(i-1)*(h-1)+(h-1))=Adiag;
end
for i=2:h-1
    A((i-2)*(h-1)+1:(i-2)*(h-1)+(h-1),(i-1)*(h-1)+1:(i-1)*(h-1)+(h-1))=Aother;
    A((i-1)*(h-1)+1:(i-1)*(h-1)+(h-1),(i-2)*(h-1)+1:(i-2)*(h-1)+(h-1))=Aother;
end
%f = -2*M*sin(M*x)*cosh(M*y);
for i=1:h-1
    for j=1:g-1
        f(i,j)=Approxfunc(i*dx,j*dy,M);
    end
end
f=reshape(f,length(f)^2,1);
u=A\(dx^2*f);
u = reshape(u,h-1,h-1);
subplot(2,1,1) 
surf(u)

for i=1:h-1
    for j=1:g-1
        uExact(i,j)=Exactfunc(i*dx,j*dy,M,Lx);
    end
end
subplot(2,1,2) 
surf(uExact)

function z = Approxfunc(x,y,M)
z=-2*M*sin(M*x)*cosh(M*y);
end
function w = Exactfunc(x,y,M,L)
w=(L-y)*sin(M*x)*sinh(M*y);
end