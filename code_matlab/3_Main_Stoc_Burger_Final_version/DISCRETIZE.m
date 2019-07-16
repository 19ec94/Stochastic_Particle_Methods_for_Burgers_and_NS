function [x,dx]= DISCRETIZE(xStart,xEnd,nb_cells)
x  = linspace(xStart,xEnd,nb_cells+1);
dx = x(2)-x(1);

end