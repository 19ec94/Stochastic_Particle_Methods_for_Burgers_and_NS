function [avg_new_cell_vel]= NEW_CELL_VEL_UPDATE (total_nb_particles,nb_cells, ...
    x,par_new,avg_new_cell_vel)
num_par_in_cell= zeros(1,nb_cells);
%position_x =0; 
for par_number=1:total_nb_particles
    %There is a scope to paralalize this loop, as it runs thorugh number of
    %particles.For each particle  it calculates the cell number and
    %accumulates in the new_cell_curl and accumulates count in num_par_in_a_cell.
    %Attention should be paid to the "new_cell_curl","num_par_in_a_cell".
    %Crtical or atomic would be fine.
    xl=1; xr=nb_cells+1; xm=floor((xl+xr)/2);
    position_x = par_new(1,par_number);
    [xl,xr] = find_cell_x(position_x,xl,xr,xm,x);
    %Note xl and yd are the actual new cell identification
    %I have to indentify the actual new_curl_cell array to store
    %the curl value
    avg_new_cell_vel(xl)=avg_new_cell_vel(xl) + par_new(2,par_number);
    %I also store the count of the number of particles in the
    %"xl+(yd-1)*nb_cells_in_x" cell. 
    num_par_in_cell(xl)= num_par_in_cell(xl)+1;
end
%New average cell curl
%I divide each cell by the current number of particles I get the average
%value of curl of the cell
avg_new_cell_vel = avg_new_cell_vel ./ num_par_in_cell;
end

