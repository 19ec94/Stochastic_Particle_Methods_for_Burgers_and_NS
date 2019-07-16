function avg_new_cell_vel = AVG_CELL_VEL_UPDATE(x,nb_cells,total_nb_particles,par_new, ...
    avg_new_cell_vel, num_par_in_a_cell)
for j=1:nb_cells
        for i=1:total_nb_particles
            if  (par_new(1,i)>= x(1,j) ) && ( par_new(1,i)< x(1,j+1) )
                avg_new_cell_vel(1,j)=avg_new_cell_vel(1,j)+par_new(2,i);
                num_par_in_a_cell = num_par_in_a_cell +1;
            end
        end
        avg_new_cell_vel(1,j)= avg_new_cell_vel(1,j)/num_par_in_a_cell;
        num_par_in_a_cell =0;
    end

end