%%Uniform random variable
a=3;
b=5;
random_vars_b_a = (b-a)*rand+a;

%%  finite Discrete Random variable generation
%% Brownian Bridge
x=100;
y=105;
T=1;
NbTraj=2;
NbPas=10;
DeltaT = T/NbPas;

dw =[zeros(NbTraj,1),sqrt(DeltaT)*randn(NbTraj,NbPas)];
w =cumsum(dw,2);
Inc =[zeros(NbTraj,1), ...
    repmat(DeltaT/T*(y-x-w(:,NbPas+1)),1,NbPas)];
rep=w+x+cumsum(Inc,2);
figure
plot([0:DeltaT:1],rep)
%% steps to compute Euler approximations
%NbTraj -number of trajectory
%Number of step to simulate 

%%
x=0:0.01:1;
u_0=sin(2*pi*x);
hold on
plot (x,u_0)
for i=0:0.05:.25
    u=sin(2*pi*x);
    x=x+i*u;
    if mod(i,0.05)==0
        hold on
        plot(x,u)
    end
end
%hold on
%plot(x,u)

%%
[r, the] = cart2plr(13,4)

