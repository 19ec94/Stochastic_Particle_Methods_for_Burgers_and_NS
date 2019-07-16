function avg_cell_vel = INIT_CELL_VEL(x,nb_cells,avg_cell_vel)
dx = x(2)-x(1);
for i=1:nb_cells
    avg_cell_vel(1,i)=sin(2*pi* (x(1,i)+(x(1,i+1)-x(1,i))/2));
    %avg_cell_vel(1,i)= exp(-2 *  ( (x(1,i)+(dx/2))^2 ) );
    
    %if (i ~= nb_cells/2)
    %    avg_cell_vel(1,i)=0;
    %else
    %    avg_cell_vel(1,i)=10;
    %end
    % ui = exp ( - 2.0 * x.^2 );
    
end
end