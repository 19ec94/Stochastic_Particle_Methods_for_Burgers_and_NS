function par_old = INIT_PAR_VEL(par_old,nb_cells,ini_nb_particle_in_a_cell,avg_cell_vel)
%This function assigns each particle a velocity. It is important to note 
%that the particles themselves don't possess any velocity. Instead it takes
%the velocity of the cell it is from. For example, initial velocity of 
%particles that are in "cell one" have velocity of "avg_cell_vel(1)"
%particles in the fifth cell will have velocity of "avg_cell_vel(5) etc.,"
for i=1:nb_cells
    for j=(((i-1)*ini_nb_particle_in_a_cell)+1) : (i*ini_nb_particle_in_a_cell)
        par_old(2,j)=avg_cell_vel(1,i);
    end
end
end