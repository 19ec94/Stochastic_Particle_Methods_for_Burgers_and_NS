%Burger's equation solution
%using Crank-Nicolson Implicit finite-difference scheme
%constants
x0 = 0;
x1 = 1;
t0 = 0;
t1 = 1;
dx = 0.05;
dt = 0.05;
alpha = 0.01;

%constants

s = alpha * dt / dx / dx;
nx = (x1-x0)/dx;
nt = (t1-t0)/dt;

%matrices.

U = zeros(nt,nx);
D = zeros(nx);
Ut = zeros(nx);
A = zeros(nt,nx);
%CONDITIONS


%initial condition.
% at t=t0, U(x,t0) =  Ut0(x)
% Ut0(x) = 100*sin(2*pi*x)

for i=1:1:nx
    x =(x0 +(i)*dx);
    U(1,i) = sin(2*pi*x);
end

%boundary conditions.
% at x=xo, U(x0,t) =  Ux0(t);
% at x=x1, U(x1,t) =  Ux1(t);
% Ux0(t) = 0 ; Ux1 = 0;

for i=1:1:nt
    %for x=x0
    U(i,1) = 0;
    %for x=x1
    U(i,nx) = 0;
end