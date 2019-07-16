clear all; 
close all;
clc;
%Variables
%nb_cells,nb_particles,dt,final_time,nt,nu,xStart,xEnd
nb_cells=100; nb_particles=1;
dt=0.01; final_time=0.1; 
nt=uint32(final_time/dt);
nu =0.0;
xStart=0; xEnd=1;
% x(nb_cells+1),dx
x = linspace(xStart,xEnd,nb_cells+1);
dx=x(2)-x(1);
%cell_centre(nb_cells),cell_vel(nb_cells),new_cell_vel(nb_cells)
cell_centre = zeros(1,nb_cells);
for i=1:nb_cells
    cell_centre(i)=x(i)+(dx/2);
end
cell_vel = 0 .* cell_centre;
new_cell_vel = 0 .* cell_vel;
for i=1:nb_cells
    cell_vel(i)=sin(2*pi*cell_centre(i));
end
%U (nt+1,nb_cells)
U = zeros(nt+1,nb_cells);
U(1,:)= cell_vel(1,:);
%total_nb_particles
total_nb_particles = nb_particles * nb_cells;
%par_old(2,total_nb_particles),par_new(3,total_nb_particles);
par_old = zeros(3,total_nb_particles);
par_new = 0 .* par_old;
%cell_number(total_nb_particles)
cell_number=0 .* par_old(1,:);

for i=1:nb_cells
    for j=(((i-1)*nb_particles)+1):(i*nb_particles)
        par_old(1,j)=cell_centre(1,i);
        par_old(2,j)=cell_vel(1,i);
        par_old(3,j)=i;
    end
end
par_new=par_old;
%%
for time = 1:nt
    for i=1:total_nb_particles
        par_new(1,i)=par_old(1,i)+(dt*par_old(2,i))+(sqrt(2*nu)* sqrt(dt)*randn);
        if(par_new(1,i) > xEnd)
            par_new(1,i)= par_new(1,i)-xEnd;
        elseif(par_new(1,i) < xStart)
            par_new(1,i)=par_new(1,i)+xEnd;
        end
    end
    for i=1:total_nb_particles
        position_x =par_new(1,i);
        xl=1;
        xr=nb_cells+1;
        xm= floor((xl+xr)/2);
        [xl,xr] = find_cell_x(position_x,xl,xr,xm,x);
        par_new(3,i)=xl;
    end
    new_cell_vel = zeros(1,nb_cells);
    particles = zeros(1,nb_cells);
    for i=1:total_nb_particles
        c = par_new(3,i) ;
        new_cell_vel(1,c) = new_cell_vel(1,c) + par_new(2,i);
        particles(1,c)= particles(1,c)+1;
    end
    for i=1:nb_cells
        if(particles(1,i)==0)
            %disp 'no particles'
            %disp(i)
            %disp(time*dt)
            %return
            new_cell_vel(1,i)=0;
        else
            new_cell_vel(1,i) = new_cell_vel(1,i) / particles(1,i);
        end
    end
    
    for i=1:total_nb_particles
        c = par_new(3,i) ;
        par_new(2,i) = new_cell_vel(1,c);
    end
    par_old = par_new;
    U(time+1,:)=new_cell_vel(1,:);
end
%%

%for i=nt+1
%    plot(cell_centre,U(i,:),'o-'); hold on
%end

%hold on
%%
grid = burgers_time_viscous_1(nb_cells,dt,final_time,nu,4,1);
%%
%for i=nt+1
%    plot(cell_centre,UN(i,:)); hold on
%end
%%
%for i=1:nt+1
%    n(i)=norm(U(i,:)-UN(i,:));
%end
%%
%T = 0:dt:final_time;
%%
%plot(T,n);

%}
%%
%
%serial=readmatrix('serial.txt');
%[m,n]=size(serial);
%serial(:,n)=[];
%%
%parallel=readmatrix('parallel.txt');
%[m,n]=size(parallel);
%parallel(:,n)=[];

%%
%for i=1:nt
    %n2(i)=norm(U(i,:)-C(i,:));
  %  if(n2(i) < (10^-8))
  %      n2(i)=0;
  %  end
%end
%norm( (parallel(nt+1,:)-serial(nt+1,:)))
%norm( (parallel(1,:)-grid(1,:)))
%norm( (serial(nt+1,:)- grid(nt+1,:)))
%norm( (parallel(1,:)-U(1,:)))
%norm( (serial(nt+1,:)-U(nt+1,:))) 
norm( (grid(nt+1,:)-U(nt+1,:))) 
%%
%plot(cell_centre,serial(nt+1,:),'o-'); hold on
%}
%plot(cell_centre,parallel(nt+1,:),'r*-'); hold on
plot(cell_centre,U(nt+1,:),'bo-');