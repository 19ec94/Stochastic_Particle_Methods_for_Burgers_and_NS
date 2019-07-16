clear all;
clc;
close all;
%Finalized version, tested against Main_Stochastic_burgers
%% Inputs
%%Time step
dt =0.001;
final_time=1.0;
ini_cond = 1;
diff_co_eff=0.01;
% We define the domain of our problem. Since it is a 1D case we only need
% a start point and an end point.
xStart=0; xEnd = 1;

% We discretize the given domain into equal number of sub-domains called
% cells. To do just that we need to mention how many cells we want
nb_cells=100;

% This function now discretize the 1D domain into "nb_cells+1" points and
% the variable "x" stores the absolute length of each point from "xStart"
% for example x(1)=xStart=0, x(2)=... , x(nb_cells+1)=xEnd=1. Additionaly
% it also returns the length of a single cell "dx".
[x,dx]= DISCRETIZE(xStart,xEnd,nb_cells);

%We define initial number of particles in each cell. We tacitly assume that
%each cell contains equal number of particles
ini_nb_particle_in_a_cell= 100;

% This function calculates total number of particles in the whole domain.
% It must be constant throughout our simulations. I.e, there is no loss of
% particles at any given time.
total_nb_particles = ini_nb_particle_in_a_cell*nb_cells;


%We create a variable "par_old" that stores position and velocity of each
%particle at time "t". We also create another exactly same variable
%"par_new" that would store the same details. we create two
%seperate variables because one will be used to store the new position and
%velocity at time "t", the other stores the old position and velocity at
% "t-1"
par_old = zeros(2,total_nb_particles);
par_new = zeros(2,total_nb_particles);

%We are interested in calculating the cell velocity. The following
%variables store the velocity of each cell at given time "t".
avg_cell_vel = zeros(1,nb_cells);
avg_new_cell_vel = zeros(1,nb_cells);

%diffusion coefficient
nt= final_time/dt;
f_u = zeros(nt+1,nb_cells);

%% We initialize the starting point for each particle
par_old = INIT_POS(x,dx,nb_cells,ini_nb_particle_in_a_cell);
%% we assign initial velocity to each cell
avg_cell_vel = INIT_CELL_VEL(x,nb_cells,avg_cell_vel,ini_cond);
f_u(1,:) = avg_cell_vel;
%% we assign initial velocity to each particle in all cells
%Very important to note that here I give the variable "par_old" to the
%funtion as an input. otherwise the program takes it as fresh value and
%overwrites the previous values by  zeros.
par_old = INIT_PAR_VEL(par_old,nb_cells,ini_nb_particle_in_a_cell,avg_cell_vel);
%% Copying the particle data from old to new
%In this crucial step we copy the position and velocity of each particles
%from "par_old(:,:)" to "par_new(:,:)". This is,in fact, very important. We may not
%need the old position because while calculating new position we overwrite
%the position value in "par_new(1,:)" but it is crucial to store the old
%velocity in the "par_new(2,:). This will tremendously simplyfy the
%calculation at later stage, in timing loop
par_new = par_old;
%% Evolution of particles according to Euler-maryama scheme
%rng('default');
tic
num_par_in_a_cell =0; t=0;
for time=1:nt
    t = t+1;
    %rng default; rng(1);
    %Update the position of particle with the increment of time step "dt"
    par_new = PAR_POS_UPDATE(par_new,par_old,dt,diff_co_eff,total_nb_particles,xStart,xEnd);
    
    %we Initialize the new average cell velocity array to zeros. If we
    %don't the values will keep on addtion time step after time step.
    avg_new_cell_vel =zeros(1,nb_cells);
    
    %Calculate new average cell velocity
    avg_new_cell_vel = AVG_CELL_VEL_UPDATE(x,nb_cells,total_nb_particles,par_new, ...
        avg_new_cell_vel, num_par_in_a_cell);
    
    %Assign the new velocity to all particles in each cell
    [par_new,par_old] = COPY_DATA(nb_cells,total_nb_particles,par_new,par_old,x,avg_new_cell_vel);
    f_u(t+1,:)=avg_new_cell_vel(1,:);
    
end
toc

%% ploting solution at Initial and Final time
if nb_cells==10000
    for i=1:nb_cells
        x_new_1(i)= x(i)+(x(i+1)-x(i))/2;
    end
    figure(1)
    set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
    plot(x_new_1,f_u(1,:))
    hold on
    plot(x_new_1,f_u(nt+1,:))
    
    % Conventional burger
    U = burgers_time_viscous_1(nb_cells,dt,final_time,diff_co_eff,4,ini_cond);
    plot(x_new_1,U(1,:))
    plot(x_new_1,U(nt+1,:))
    
    %% Surface plotting for all time-Stocha
    T=0:dt:final_time;
    figure(2)
    set(gcf,'units','normalized','position',[0.25 0.52 0.23 0.32]);
    surf(x_new_1,T,f_u)
    shading interp
    xlabel('x'), ylabel('t'), zlabel ('Stocha u(x,t)');
    grid on
    %colormap('Gray');
    %% Surface plotting for all time - conven
    T=0:dt:final_time;
    figure(3)
    set(gcf,'units','normalized','position',[0.49 0.52 0.23 0.32]);
    surf(x_new_1,T,U)
    shading interp
    xlabel('x'), ylabel('t'), zlabel ('Conventional U(x,t)');
    grid on
    %colormap('Gray');
    %% Norm Plotting  for all time
    for i=1:(nt+1)
        n(i)=norm(f_u(i,:)-U(i,:));
    end
    figure(4)
    set(gcf,'units','normalized','position',[0.73 0.52 0.23 0.32]);
    plot(T,n);
end

%% Plotting for case of 1000 grid cells
if nb_cells==10000
     for i=1:nb_cells
        x_new_1(i)= x(i)+(x(i+1)-x(i))/2;
    end
    for i=1:(nb_cells/100)
        x_new_100(i) = x_new_1(i*100);
    end
    for i=1:(nt+1)
        for j=1:(nb_cells/100)
            f_u_100(i,j)= f_u(i,j*100);
        end
    end
    
    figure(1)
    set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
    plot(x_new_100,f_u_100(1,:))
    hold on
    plot(x_new_100,f_u_100(nt+1,:))
    
    % Conventional burger
    U = burgers_time_viscous_1(100,dt,final_time,diff_co_eff,4,ini_cond);
    plot(x_new_100,U(1,:))
    plot(x_new_100,U(nt+1,:))
    
    %% Surface plotting for all time-Stocha
    T=0:dt:final_time;
    figure(2)
    set(gcf,'units','normalized','position',[0.25 0.52 0.23 0.32]);
    surf(x_new_100,T,f_u_100)
    shading interp
    xlabel('x'), ylabel('t'), zlabel ('Stocha u(x,t)');
    grid on
    %colormap('Gray');
    %% Surface plotting for all time - conven
    T=0:dt:final_time;
    figure(3)
    set(gcf,'units','normalized','position',[0.49 0.52 0.23 0.32]);
    surf(x_new_100,T,U)
    shading interp
    xlabel('x'), ylabel('t'), zlabel ('Conventional U(x,t)');
    grid on
    %colormap('Gray');
    %% Norm Plotting  for all time
    for i=1:(nt+1)
        n(i)=norm(f_u_100(i,:)-U(i,:));
    end
    figure(4)
    set(gcf,'units','normalized','position',[0.73 0.52 0.23 0.32]);
    plot(T,n);
end
for i=1:nb_cells
        x_new_1(i)= x(i)+(x(i+1)-x(i))/2;
    end
    figure(1)
    set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
    plot(x_new_1,f_u(1,:))
    hold on
    plot(x_new_1,f_u(nt+1,:))
      U = burgers_time_viscous_1(nb_cells,dt,final_time,diff_co_eff,4,ini_cond);
    plot(x_new_1,U(1,:))
    plot(x_new_1,U(nt+1,:))