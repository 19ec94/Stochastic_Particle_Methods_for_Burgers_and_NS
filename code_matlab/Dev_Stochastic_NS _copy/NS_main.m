clc;
clear;
%% Domain
dt =0.001; final_time=0.5;
%We want to solve NS equation in a 2D domain. For that we define our domain.
%we specify "x and y" domain.
%Start and End point in x and y direction
xStart = 0; xEnd = 1;  yStart = 0; yEnd= 1;
%We want to discretize the "x" and "y" domain into "N" and "N" number of
% sub-domians. Importantly we can have different number of cells in "x" and
%"y" directions. So we need two variables.
nb_cells_in_x = 100 ; nb_cells_in_y = 100;
[total_nb_cells,x_domain,y_domain,dx,dy]=DISCRETIZE (xStart, ...
    xEnd,yStart,yEnd,nb_cells_in_x,nb_cells_in_y);
%We store each cell's absolute length in "x" and "y" directions from the
%respective origin
[x_cell,y_cell] = CELL_LENGTH(nb_cells_in_x,nb_cells_in_y,x_domain, ...
    y_domain);
%Number of particle in a cell
nb_of_particles_in_a_cell = 100;
%Total_number_of_particles
total_nb_particles= nb_of_particles_in_a_cell*nb_cells_in_x*nb_cells_in_y;
%populate the cells with particles. It is  5*number_of_particles array
%that stores the whole data about the particle's  x,y,
%velocity_x,velocity_y,vorticity  of the cell it comes from
%par_old = zeros(5,total_nb_particles);
par_new = zeros(5,total_nb_particles);
%diffusion coefficient
diff_co_eff=0.005;
%Each cell has a 2D velocity vector. And particles have the velocity of the
%cell it comes from.  We define two variables to store the values at time
%"t-1" and "t"
cell_vel = zeros(2,total_nb_cells);
new_cell_vel = zeros(2,total_nb_cells);
%cell_curl = zeros(1,total_nb_cells);
new_cell_curl = zeros(1,total_nb_cells);
current_cell = zeros(1,total_nb_particles);
%time_step
num_par_in_a_cell=0;
%%
%We store each cells xStart,xEnd, yStart and yEnd and it's cell number.
%We can use it later ?
%cell_coord = zeros(5,total_nb_cells);
%We store the centre point of each cell. This is used in initializing
%particles at this position and calculating velocity and curl at this point.
%cell_centre_coord = zeros(2,total_nb_cells);
[cell_coord, cell_centre_coord] = ...
    CELL_COORDINATES(x_cell,x_domain,y_cell,y_domain,nb_cells_in_x,...
    nb_cells_in_y,total_nb_cells);
%% Initialize each particle's position in the domain
%We set each particle's initial position as centre of the respective cell
par_old = INIT_POS_PAR(total_nb_cells, nb_of_particles_in_a_cell, ...
    nb_cells_in_x,dx,dy,cell_centre_coord);
%% initial velocity
%we calculate cell velocity at it's centre position
for i=1:total_nb_cells
    cell_vel(1,i)=sin(2*pi* cell_centre_coord(2,i));
    cell_vel(2,i)=cos(2 *pi*cell_centre_coord(1,i));
end
%%  Initial cell curls
%we calculate cell curl at it's centre position
cell_curl = CURL_FUNCTION(total_nb_cells,cell_centre_coord);
%% Initialize each particle's velocity and vorticity
for i = 1:total_nb_cells
    for xInd=(((i-1)*nb_of_particles_in_a_cell)+1):(i*nb_of_particles_in_a_cell)
        par_old(3,xInd)=cell_vel(1,i);
        par_old(4,xInd)=cell_vel(2,i);
        par_old(5,xInd)=cell_curl(1,i);
    end
end
%% copy
par_new(1,:) = par_old(1,:);
par_new(2,:) = par_old(2,:);
par_new(5,:) = par_old(5,:);
t=0;
dist = zeros(2,total_nb_cells);
pt = zeros(2,total_nb_cells);
%%
for time=0:dt:final_time
    t=t+1;
    % update the particle position
    [par_new] = POS_UPDATE(total_nb_particles,par_old,par_new,diff_co_eff,dt,...
        xStart,xEnd,yStart,yEnd);
    
    % update cell vorticity
    [current_cell,new_cell_curl]= NEW_CELL_CURL (current_cell,total_nb_particles,total_nb_cells,nb_cells_in_x, ...
        nb_cells_in_y,x_domain,y_domain,par_new,new_cell_curl);
    
    % Update cell velocity
    %cell velocity after calculating contribution from all the cells
    %{
    for i=1:total_nb_cells
        pt(1,i)=cell_centre_coord(1,i)+dx;
        pt(2,i)=cell_centre_coord(2,i)+dy;
        for j=1:total_nb_cells
            dist(1,j) = pt(1,i) - cell_centre_coord(1,j);
            dist(2,j) = pt(2,i) - cell_centre_coord(2,j);
            modulus_dist = sqrt(dist(1,j)^2+dist(2,j)^2);
            new_cell_vel(1,i)= new_cell_vel(1,i) + ((new_cell_curl(1,j)*dist(1,j))./modulus_dist);
            new_cell_vel(2,i)= new_cell_vel(2,i) + ((new_cell_curl(1,j)*dist(2,j))./modulus_dist);
        end
        new_cell_vel(1,i) = (1/2*pi)*new_cell_vel(1,i);
        new_cell_vel(2,i) = (1/2*pi)*new_cell_vel(2,i);
    end
    %}
    
    %cell velocity after celculating only from the respective curl
    new_cell_vel = 0 .* new_cell_vel;
    for i=1:total_nb_cells
        pt(1,i)=cell_centre_coord(1,i)+(dx/10);
        pt(2,i)=cell_centre_coord(2,i)+(dy/10);
        dist(1,i) = cell_centre_coord(1,i)-pt(1,i);
        dist(2,i) = cell_centre_coord(2,i)-pt(2,i);
        modulus_dist = sqrt(dist(1,i).^2+dist(2,i).^2);
        new_cell_vel(1,i)= new_cell_vel(1,i) + ((new_cell_curl(1,i)*dist(1,i))./modulus_dist);
        new_cell_vel(2,i)= new_cell_vel(2,i) + ((new_cell_curl(1,i)*dist(2,i))./modulus_dist);
        
        new_cell_vel(1,i) = (1/2*pi)*new_cell_vel(1,i);
        new_cell_vel(2,i) = (1/2*pi)*new_cell_vel(2,i);
    end
    %Update newposition to old position and new average velocity to old par
    %and new par
    for i=1:total_nb_particles
        r = current_cell(1,i);
        par_old(3,i)= new_cell_vel(1,r);
        par_old(4,i)= new_cell_vel(2,r);
        par_new(5,i)= new_cell_curl(1,r);
    end
end

%%
%%
x = cell_centre_coord(1,:);
y = cell_centre_coord(2,:); 
%%
u = new_cell_vel(1,:);
v=new_cell_vel(2,:);
quiver(x,y,u,v);
%%
u1 = cell_vel(1,:);
v1=cell_vel(2,:);
quiver(x,y,u1,v1);
