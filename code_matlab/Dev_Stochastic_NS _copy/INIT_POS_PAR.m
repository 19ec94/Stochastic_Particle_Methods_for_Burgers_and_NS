function par_old = INIT_POS_PAR(total_nb_cells, nb_of_particles_in_a_cell, ...
    nb_cells_in_x,dx,dy,cell_centre_coord)


for i=1:total_nb_cells
    for j=(((i-1)*nb_of_particles_in_a_cell)+1):(i*nb_of_particles_in_a_cell)
        par_old(1,j)=cell_centre_coord(1,i);
    end
end
for i=1:total_nb_cells
    for j=(((i-1)*nb_of_particles_in_a_cell)+1):(i*nb_of_particles_in_a_cell)
        par_old(2,j)=cell_centre_coord(2,i);
    end
end
end