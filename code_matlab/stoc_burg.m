%%Domain x=[0,1]
x =0:1/100:1;
%Choose a point in the domain
y=x(ceil(end/2));
%choose number of particles
NbTraj =10;
%Number of time steps for each particles ?????
NbPas = 500;
%Create an array for all the particles for all the time steps
w=zeros(NbTraj,NbPas);
%Assign the intial value for all the particle
w(:,1) =y;
ini_u = sin(2*pi*y);
%Time step
dt = 0.01;
%Let the particle evolve one by one
for particle=1:1:NbTraj
    for i=1:1:NbPas
        w(particle,i+1)=w(particle,i)+ dt * sin(2*pi*w(particle,i)) + sqrt(2*0.1)* sqrt(dt)* randn;
        if (w(particle,i+1)>=1)  || ( w(particle,i+1) <= 0 )
            break;
        end
    end
end
final_u= sin(2*pi*mean(w(:,20)));
%x = linspace(0,1,length(w));
%plot(x,w(:,1),'*')
%hold on
%plot(x,w(:,20),'*-')
