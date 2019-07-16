clc;
clear;
%% Domain
%Start and ENd point in x direction
xStart = 0; xEnd = 1;
%Start and End point in y direction
yStart = 0; yEnd= 1;
%number of intervals(Cells) in x and y direction
nb_cells_in_x = 10; nb_cells_in_y = 10;
%total number of cells
total_nb_cells = nb_cells_in_x* nb_cells_in_y;
%discretize the domain into nb_cells+1 points
x_domain= linspace(xStart,xEnd,nb_cells_in_x+1);
y_domain= linspace(yStart,yEnd,nb_cells_in_y+1);
dx = x_domain(2)-x_domain(1);
dy = y_domain(2)-y_domain(1);
%x_cell stores the absolute length in x direction from xstart
x_cell = zeros(1,nb_cells_in_x);
%y_cell stores the absolute length in y direction from ystart
y_cell = zeros(1,nb_cells_in_y);
for index=1:nb_cells_in_x
    x_cell(index)=x_domain(index+1);
end
for index=1:nb_cells_in_y
    y_cell(index)=y_domain(index+1);
end
%Number of particle in a cell
nb_of_particles_in_a_cell = 1;
%Total_number_of_particles
nb_particles= nb_of_particles_in_a_cell*nb_cells_in_x*nb_cells_in_y;
%populate the cells with particles. It is  4*number_of_particles array
%that stores the whole data about the particle. x,y,,cell number, velocity
%of the cell it come from
par_old = zeros(3,nb_particles);
par_new = zeros(3,nb_particles);
%diffusion coefficient
diff_co_eff=0.005;
cell_vel = zeros(2,total_nb_cells);
new_cell_vell = zeros(2,total_nb_cells);
%%
cell_coord = zeros(2,total_nb_cells);
cell_centre_coord = zeros(2,total_nb_cells);
for i=1:nb_cells_in_x %pointer is at individual column and progress through row
    for j=i:nb_cells_in_x:total_nb_cells
        cell_coord(1,j)=x_cell(1,i);
        cell_centre_coord(1,j)=x_domain(1,i)+ (x_domain(1,i+1)-x_domain(1,i))/2;
    end
end
for i=1:nb_cells_in_y %pointer is at row and progress through row
    for j= (((i-1)*nb_cells_in_x)+1):(i*nb_cells_in_x)
        cell_coord(2,j)=y_cell(1,i);
        cell_centre_coord(2,j)=y_domain(1,i)+ (y_domain(1,i+1)-y_domain(1,i))/2;
    end
end
%% Initialize each particle's position in the domain
for i=1:total_nb_cells
    for j=(((i-1)*nb_of_particles_in_a_cell)+1):(i*nb_of_particles_in_a_cell)
        if (mod(i,nb_cells_in_x) ~=1)
            par_old(1,j) = cell_coord(1,i-1)+ (dx/2);
        else
            par_old(1,j) = cell_coord(1,i)/2 ;
        end
    end
end
for i=1:total_nb_cells
    for j=(((i-1)*nb_of_particles_in_a_cell)+1):(i*nb_of_particles_in_a_cell)
        if (i> nb_cells_in_x)
            par_old(2,j) = cell_coord(2,i-nb_cells_in_x)+ (dy/2);
        else
            par_old(2,j) = cell_coord(2,i)/2 ;
        end
    end
end

%% initial velocity
for i=1:total_nb_cells
    cell_vel(1,i)=sin(cell_centre_coord(2,i));
    cell_vel(2,i)=cos(cell_centre_coord(1,i));
end

%% cell vorticity

mesh_x = zeros(1,total_nb_cells);
mesh_y = zeros(1,total_nb_cells);
new_curl=zeros(1,nb_cells_in_x*nb_cells_in_y);
%u = zeros(1,total_nb_cells);
%v = zeros(1,total_nb_cells);
%dummy_curl =zeros(1,total_nb_cells);
for i=1:total_nb_cells
    mesh_x(i) = cell_centre_coord(1,i);
    mesh_y(i) = cell_centre_coord(2,i);
end
[X,Y]=meshgrid(mesh_x,mesh_y);
%[X,Y] = meshgrid(x_domain,y_domain);
%Vxx = Y.*X.^2+3.*Y.^2; Vyy = 2.*X.*Y+X.^2; CURL=-X.^2+2.*X-4.*Y;
%u = sin(Y); v=cos(X); 
CURL= -sin(X)-cos(Y);
%[dummy_curl]=curl(X,Y,u,v);
%subplot(2,1,1),contourf(X,Y,CURL),colorbar
%subplot(2,1,2),contourf(X,Y,dummy_curl),colorbar
for i=1:(nb_cells_in_x*nb_cells_in_y)
    new_curl(i) = CURL(i,i);
end
%% 
curl_third = zeros(1,total_nb_cells);
ori_curl_third = zeros(1,total_nb_cells);
syms f(x,y,z)
f(x,y,z) = [sin(y) cos(x) 0];
sm_curl = curl(f,[x y z]);
func_sm_curl = matlabFunction(sm_curl,'vars', [x y z],'Optimize',false);
for i=1:total_nb_cells
    numerical_sm_curl = func_sm_curl(cell_centre_coord(1,i), ...
         cell_centre_coord(2,i),0); 
    curl_third(i) = numerical_sm_curl(3);
    ori_curl_third(i) = -sin(cell_centre_coord(1,i))- ...
        cos(cell_centre_coord(2,i));
end
%%





