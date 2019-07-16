function [cell_coordinates, cell_centre_coordinates] = ...
    coordinates(x_cell,x_domain,y_cell,y_domain,nb_cells_in_x,...
    nb_cells_in_y,total_nb_cells)
for i=1:nb_cells_in_x %pointer is at individual column and progress through row
    for j=i:nb_cells_in_x:total_nb_cells
        cell_coordinates(1,j)=x_cell(1,i);
        cell_centre_coordinates(1,j)=x_domain(1,i)+ (x_domain(1,i+1)-x_domain(1,i))/2;
    end
end
for i=1:nb_cells_in_y %pointer is at row and progress through row
    for j= (((i-1)*nb_cells_in_x)+1):(i*nb_cells_in_x)
        cell_coordinates(2,j)=y_cell(1,i);
        cell_centre_coordinates(2,j)=y_domain(1,i)+ (y_domain(1,i+1)-y_domain(1,i))/2;
    end
end


end
