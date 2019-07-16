function ui = ic_gaussian ( x,ini_cond )

%*****************************************************************************80
%
%% IC_GAUSSIAN evaluates the initial condition for a Gaussian.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    22 April 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X(*), the node coordinates.
%
%    Output, real UI(*), the value of the initial condition at each node.
%
%ui = exp ( - 2.0 * x.^2 );
dx=x(2)-x(1);
nb_cells = length(x)-1;
if ini_cond==1
    for i=1:nb_cells
        ui(1,i)=sin(2*pi* (x(1,i)+(dx/2)) );
    end
elseif ini_cond==2
    for i=1:nb_cells
        ui(1,i)= exp(-2 *  ( (x(1,i)+(dx/2))^2 ) );
    end
elseif ini_cond==3
    for i=1:nb_cells
        if (i ~= nb_cells/2)
            ui(1,i)=0;
        else
            ui(1,i)=10;
        end
    end
else
    disp 'initial condition failed'
end
return
end
