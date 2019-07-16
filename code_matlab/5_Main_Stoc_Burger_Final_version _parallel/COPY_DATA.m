function [par_new,par_old] = COPY_DATA(nb_cells,total_nb_particles,par_new,par_old,x,avg_new_cell_vel)
for j=1:nb_cells
    for i=1:total_nb_particles
        if  (par_new(1,i)>= x(1,j) ) && ( par_new(1,i)<= x(1,j+1) )
            par_old(1,i) = par_new(1,i);
            par_old(2,i) = avg_new_cell_vel(1,j);
            par_new(2,i) = avg_new_cell_vel(1,j);
        end
    end
end

end