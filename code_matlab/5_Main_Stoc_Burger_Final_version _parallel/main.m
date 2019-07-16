
%Finalized version, tested against Main_Stochastic_burgers
%% Inputs
%%Time step
dt =0.0001;
final_time=1.0;
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
diff_co_eff=0.005;



%% We initialize the starting point for each particle
par_old = INIT_POS(x,dx,nb_cells,ini_nb_particle_in_a_cell);
%% we assign initial velocity to each cell
avg_cell_vel = INIT_CELL_VEL(x,nb_cells,avg_cell_vel);
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
num_par_in_a_cell =0;
for time=0:dt:final_time
    %rng default; rng(1);
    %Update the position of particle with the increment of time step "dt"
    par_new = PAR_POS_UPDATE(par_new,par_old,dt,diff_co_eff,total_nb_particles,xStart,xEnd);
    
    %we Initialize the new average cell velocity array to zeros. If we
    %don't the values will keep on addtion time step after time step.
    avg_new_cell_vel =zeros(1,nb_cells);
    
    %Calculate new average cell velocity
    %avg_new_cell_vel = AVG_CELL_VEL_UPDATE(x,nb_cells,total_nb_particles,par_new, ...
    %   avg_new_cell_vel, num_par_in_a_cell);
    [avg_new_cell_vel]= NEW_CELL_VEL_UPDATE (total_nb_particles,nb_cells, ...
    x,par_new,avg_new_cell_vel);
    
    %Assign the new velocity to all particles in each cell
    [par_new,par_old] = COPY_DATA(nb_cells,total_nb_particles,par_new,par_old,x,avg_new_cell_vel); 
    
end
toc
%% ploting
for i=1:nb_cells
    x_new_1(i)= x(i)+(x(i+1)-x(i))/2;
end
plot(x_new_1,avg_cell_vel); hold on
plot(x_new_1,avg_new_cell_vel);


 