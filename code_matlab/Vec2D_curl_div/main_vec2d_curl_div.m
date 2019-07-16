%%
syms x y z real
F = [ cos(x+2*y), sin(x-2*y) ];
vectorfield(F,-2:.25:2,-2:.25:2)
hold on

%%
G = curl([F,0],[x y z])                % need vector field of 3 components for curl

ezpcolor(G(3),[-2.5 2.5 -2.5 2.5]); hold on
vectorfield(F,-2:.25:2,-2:.25:2); hold off
colorbar; 
colorposneg
title('colors show curl: blue is clockwise, red is counterclockwise')
%%
clc; clear;
x = [0.1 0.2];
y = [0.1 0.2];
[X,Y]=meshgrid(x);
Vxx = X; Vyy=Y.^2;
[curlz]=curl(X,Y,Vyy,Vxx);
%%
clc;clear;
N=16;
[X,Y]=meshgrid(linspace(-4,4,N),linspace(-4,4,N)); 
u =Y.*X.^2+3*Y.^2;
v =2*X.*Y+X.^2;
% analitic curl
CURL=-X.^2+2*X-4*Y;
subplot(2,1,1),contourf(X,Y,CURL),colorbar
% using  matlab function (curl)
[cur,va]=curl(X,Y,u,v);
subplot(2,1,2),contourf(X,Y,cur),colorbar

%%
clc;clear;
N=100;
[X,Y]=meshgrid(linspace(0,1,N),linspace(0,1,N)); 
u = sin(Y);
v =cos(X);
% analitic curl
CURL=-sin(X)-cos(Y);
subplot(2,1,1),contourf(X,Y,CURL),colorbar
% using  matlab function (curl)
[cur,va]=curl(X,Y,u,v);
subplot(2,1,2),contourf(X,Y,cur),colorbar

%%
g = divergence(F,[x y])                % find divergence

ezpcolor(g,[-2.5 2.5 -2.5 2.5]); hold on
vectorfield(F,-2:.25:2,-2:.25:2); hold off
colorbar; 
colorposneg 
title('colors show divergence: blue is sink, red is source')