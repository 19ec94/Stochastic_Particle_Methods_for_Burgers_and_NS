function avg_cell_vel = INIT_CELL_VEL(x,nb_cells,avg_cell_vel,ini_cond)
dx = x(2)-x(1);
if ini_cond==1
    for i=1:nb_cells
        avg_cell_vel(1,i)=sin(2*pi* (x(1,i)+(dx/2)));
    end
elseif ini_cond==2
    for i=1:nb_cells
        avg_cell_vel(1,i)= exp(-2 *  ( (x(1,i)+(dx/2))^2 ) );
    end
elseif ini_cond==3
    for i=1:nb_cells
        if (i ~= nb_cells/2)
            avg_cell_vel(1,i)=0;
        else
            avg_cell_vel(1,i)=10;
        end
    end
else
    disp 'initial condition failed'
end
end