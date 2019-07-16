function par_old = INIT_POS(x,dx,nb_cells,ini_nb_particle_in_a_cell)
%This Function initialize the position of all particles in each cell at the
%cell's mid-point
for i=1: nb_cells
    for j=(((i-1)*ini_nb_particle_in_a_cell)+1) : (i*ini_nb_particle_in_a_cell)
        par_old(1,j)= x(1,i)+(x(1,i+1)-x(1,i))/2;
    end
end
