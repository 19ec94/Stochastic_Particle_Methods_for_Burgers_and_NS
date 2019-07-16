%clear;
%clc;
%% MATLAB code for Solution algortihm for stochastics burgers equation
%Define a domain and discretie the domain
%start point of domain
xStart=0;
%End point of domain
xEnd = 1;
%Number of cells in the domain
nb_cells=100;
%divide the domains into "cell+1" points
x = linspace(xStart,xEnd,nb_cells+1);
dx=x(2)-x(1);
%allocate space for each cell-it stores the length of each cell
cell =zeros(1,nb_cells);
for index=1:nb_cells
    cell(index)= x(index+1);
end
%Initial number of particles in each cell
particle_in_a_cell= 1000;
% total number of particles
nb_particles = particle_in_a_cell*numel(cell);
%populate domain with particles
par_old = zeros(3,nb_particles);
par_new = zeros(3,nb_particles);
%par_velocity_old = zeros(1,numel(cell));
cell_vel = zeros(1,numel(cell));
new_cell_vel = zeros(1,numel(cell));
%diffusion coefficient
diff_co_eff=0.01;
%% Initialize the starting point for each particle
%Initially all particles in the respective cells starts at the mid-point of
%the cell
for i=1:numel(cell)
    for j=(((i-1)*particle_in_a_cell)+1) : (i*particle_in_a_cell)
        par_old(1,j)= x(1,i)+(x(1,i+1)-x(1,i))/2;
        par_old(2,j)=i;
    end
end
%% Assign initial velocity to each particle
%we assign initial velocity to each particle
%for i=1:nb_particles
   % par_old(3,i)=sin(2*pi*par_old(1,i));
%end
%for i=1:numel(cell)
%    for j=(((i-1)*particle_in_a_cell)+1) : (i*particle_in_a_cell)
%        if i ~= (numel(cell)/2)
%            par_old(3,j)=0;
%        else
%            par_old(3,j)=10;
%        end
%    end
%end
%for i=1:nb_particles
%    par_old(3,i) = exp(-2*(par_old(1,i) - 1).^2);
%end

%% Calculate cell inital average velocity of each cell
for i=1:numel(cell)
    for j=(((i-1)*particle_in_a_cell)+1) : (i*particle_in_a_cell)
        cell_vel(1,i)= cell_vel(1,i) + par_old(3,j);
    end
    cell_vel(1,i) = cell_vel(1,i)/particle_in_a_cell;
end

%%
par_new = par_old;
%% Evolution of particles in one time step according to Euler-maryama scheme
dt =0.001;
velocity=0;
num_par_in_a_cell =0;
for time=0:dt:0.15
    for i=1:nb_particles
        %for j=1:nb_cells
        %    if  (par_old(1,i)>= x(1,j) ) && ( par_old(1,i)<= x(1,j+1) )
        %        velocity = par_old(3,i);
        %    end
        %end
        par_new(1,i)=par_old(1,i)+dt * par_old(3,i) + sqrt(2*diff_co_eff)* sqrt(dt)*randn;
        %periodic boundary condition
        if par_new(1,i) > xEnd
            par_new(1,i) = par_new(1,i)-xEnd;
        elseif par_new(1,i) < xStart
            par_new(1,i) = par_new(1,i)+xEnd;
        end
    end
    %%
    %sorting of new position
    % par_new_sort = ((sortrows(par_new',1)))';
    new_cell_vel =zeros(1,numel(cell));
    %%
    %%Calculate new average cell velocity
    for j=1:nb_cells
        for i=1:nb_particles
            if  (par_new(1,i)>= x(1,j) ) && ( par_new(1,i)< x(1,j+1) )
                new_cell_vel(j)=new_cell_vel(j)+par_new(3,i);
                num_par_in_a_cell = num_par_in_a_cell +1;
            end
        end
        new_cell_vel(j)= new_cell_vel(j)/num_par_in_a_cell;
        num_par_in_a_cell =0;
    end
    %%
    %%Assign the new velocity to all particles in each cell
    for j=1:nb_cells
        for i=1:nb_particles
            if  (par_new(1,i)>= x(1,j) ) && ( par_new(1,i)<= x(1,j+1) )
                par_old(1,i) = par_new(1,i);
                par_old(2,i)=j;
                par_new(3,i)=new_cell_vel(j);
            end
        end
    end
end
%% ploting
y = linspace(0,1,length(cell_vel));
plot(y,cell_vel); hold on
plot(y,new_cell_vel,'o-');
%legend('stoc-final','stoc-initial');
