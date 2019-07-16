
%%fvm
clear all
clc

  n = 20;

for i = 0:1/n:2

    u_0 = sin(2*pi*i);
    t = 0:0.05:0.5;
    x_t = i + u_0*t;
    hold on
    plot(x_t,t,'*-')
end

%% using stochastics 

y=0.5
u0=sin(2*pi*y);

