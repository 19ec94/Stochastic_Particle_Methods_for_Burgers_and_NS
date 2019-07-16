nu=0.01;
xn=10;
x=0:1/xn:1;
tn=10;
t=0:1/tn:1;

u = zeros ( xn, tn );

  for j = 1 : tn

    for i = 1 : xn

      a = ( x(i) - 4.0 * t(j) );
      b = ( x(i) - 4.0 * t(j) - 2.0 * pi );
      c = 4.0 * nu * ( t(j) + 1.0 );
      phi = exp ( - a ^ 2 / c ) + exp ( - b ^ 2 / c );
      dphi = - 2.0 * a * exp ( - a ^ 2 / c ) / c ...
             - 2.0 * b * exp ( - b ^ 2 / c ) / c;
      u(i,j) = 4.0 - 2.0 * nu * dphi / phi;
      

    end

  end
%%

