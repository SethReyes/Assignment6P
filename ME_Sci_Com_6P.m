%Seth Reyes
clear all, clc,close all

Lx= pi ;
Ly= pi ;
u0= 0 ;
uL= 0 ;
v0= 0 ;
vL= 0 ;
h=50; g=h;
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
for i=1:h-1
    for j=1:g-1
        f(i,j)=Approxfunc(i*dx,j*dy,M);
    end
end
f=reshape(f,length(f)^2,1);
u=A\(dx^2*f);
u = reshape(u,h-1,h-1);
figure('Name','Numerical (Top) vs Exact (Bottom)')
subplot(2,1,1)
title('Numerical Solution')
surf(u)
xlabel('x')
ylabel('y')
zlabel('u')

for i=1:h-1
    for j=1:g-1
        uExact(i,j)=Exactfunc(i*dx,j*dy,M,Lx);
    end
end
subplot(2,1,2)
title('Exact Solution')
surf(uExact)
xlabel('x')
ylabel('y')
zlabel('u')

z1=columnPull(Ly/4,u,uExact,Ly,g);
y=linspace(0,Ly,g+1);
y= y(2:g);
figure(2)
subplot(2,1,1)
plot(y,z1(:,1),y,z1(:,2))
legend('Numerical','Theoretical')
xlabel('x')
ylabel('u')
title('Ly/4')
z2=columnPull(Ly/2,u,uExact,Ly,g);
subplot(2,1,2)
plot(y,z2(:,1),y,z2(:,2))
legend('Numerical','Theoretical')
xlabel('x')
ylabel('u')
title('Ly/2')
z3=rowPull(3*Lx/4,u,uExact,Lx,h);
x=linspace(0,Lx,h+1);
x= x(2:h);
figure(3)
subplot(2,1,1)
plot(x,z3(:,1),x,z3(:,2))
legend('Numerical','Theoretical')
xlabel('y')
ylabel('u')
title('3Lx/4')
subplot(2,1,2)
z4=rowPull(2*Lx/4,u,uExact,Lx,h);
plot(y,z4(:,1),y,z4(:,2))
legend('Numerical','Theoretical')
xlabel('y')
ylabel('u')
title('Lx/2')
%%Functions
function z1 = columnPull(col,u,uExact,Ly,g)
z1=[u(:,round(col/Ly*(g-1))),uExact(:,round(col/Ly*(g-1)))];
end
function z2 = rowPull(col,u,uExact,Lx,h)
z2=[u(round(col/Lx*(h-1)),:),uExact(round(col/Lx*(h-1)),:)];
z2=reshape(z2,h-1,2);
end
function z = Approxfunc(x,y,M)
z=-2*M*sin(M*x)*cosh(M*y);
end
function w = Exactfunc(x,y,M,L)
w=(L-y)*sin(M*x)*sinh(M*y);
end