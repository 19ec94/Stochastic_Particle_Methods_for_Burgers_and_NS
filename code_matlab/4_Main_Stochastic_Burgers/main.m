%clear;
%clc;
%Finalized version, tested against Main_stoc_burgers_final_version
%% MATLAB code for Solution algortihm for stochastics burgers equation
final_time=1.0;
dt =0.001;
%Define a domain and discretie the domain
%start point of domain
xStart=0;
%End point of domain
xEnd = 1;
%Number of cells in the domain
nb_cells=100;
%divide the dnomains into "cell+1" points
x = linspace(xStart,xEnd,nb_cells+1);
dx=x(2)-x(1);
%allocate space for each cell-it stores the length of each cell
cell =zeros(1,nb_cells);
for index=1:nb_cells
    cell(index)= x(index+1);
end
%Initial number of particles in each cell
particle_in_a_cell= 100;
% total number of particles
nb_particles = particle_in_a_cell*nb_cells;
%populate domain with particles
par_old = zeros(2,nb_particles);
par_new = zeros(2,nb_particles);
%par_velocity_old = zeros(1,numel(cell));
cell_vel = zeros(1,nb_cells);
new_cell_vel = zeros(1,nb_cells);
%diffusion coefficient
diff_co_eff=0.01;

%% Initialize the starting point for each particle
%Initially all particles in the respective cells starts at the mid-point of
%the cell
for i=1:numel(cell)
    for j=(((i-1)*particle_in_a_cell)+1) : (i*particle_in_a_cell)
        par_old(1,j)= x(1,i)+(x(1,i+1)-x(1,i))/2;
    end
end

%%
%we assign initial velocity to each cell
for i=1:nb_cells
    cell_vel(1,i)=sin(2*pi* (x(1,i)+(x(1,i+1)-x(1,i))/2));
    %cell_vel(1,i)= exp(-2*((x(1,i)+(x(1,i+1)-x(1,i))/2) - 1).^2);
    %if (i ~= nb_cells/2)
    %    cell_vel(1,i)=0;
    %else
    %    cell_vel(1,i)=10;
    %end
end

for i=1:numel(cell)
    for j=(((i-1)*particle_in_a_cell)+1) : (i*particle_in_a_cell)
        par_old(2,j)=cell_vel(1,i);
    end
end
%%
par_new = par_old;
%% Evolution of particles in one time step according to Euler-maryama scheme

velocity=0;
num_par_in_a_cell =0;
tic
for time=0:dt:final_time
    %rng default; rng(1);
    for i=1:nb_particles
        par_new(1,i)=par_old(1,i)+dt * par_old(2,i) + sqrt(2*diff_co_eff)* sqrt(dt)*randn;
        %periodic boundary condition
        if par_new(1,i) > xEnd
            par_new(1,i) = par_new(1,i)-xEnd;
        elseif par_new(1,i) < xStart
            par_new(1,i) = par_new(1,i)+xEnd;
        end
    end
    %% sorting of new position
    new_cell_vel =zeros(1,nb_cells);
    %% Calculate new average cell velocity
    for j=1:nb_cells
        for i=1:nb_particles
            if  (par_new(1,i)>= x(1,j) ) && ( par_new(1,i)< x(1,j+1) )
                new_cell_vel(1,j)=new_cell_vel(1,j)+par_new(2,i);
                num_par_in_a_cell = num_par_in_a_cell +1;
            end
        end
        new_cell_vel(1,j)= new_cell_vel(1,j)/num_par_in_a_cell;
        num_par_in_a_cell =0;
    end
    %% Assign the new velocity to all particles in each cell
    for j=1:nb_cells
        for i=1:nb_particles
            if  (par_new(1,i)>= x(1,j) ) && ( par_new(1,i)<= x(1,j+1) )
                par_old(1,i) = par_new(1,i);
                par_old(2,i) = new_cell_vel(1,j);
                par_new(2,i) = new_cell_vel(1,j);
            end
        end
    end
end
toc
%% ploting
for i=1:nb_cells
    x_new_1(i)= x(i)+(x(i+1)-x(i))/2;
end
plot(x_new_1,cell_vel); hold on
plot(x_new_1,new_cell_vel);