function ui = uinit( x, ictype,xend )
xshift = 0;
if (ictype==1) %Shock (uL > uR)
    uL = 1;
    uR = 0;
    ui = uR + (uL-uR) * ((x-xshift) <= 0.0);
elseif (ictype==2) %Expansion (uL < uR)
    uL = 0.5;
    uR = 1.0;
    ui = uR + (uL-uR) * ((x-xshift) <= 0.0);
elseif (ictype==3) %Gaussian
    ui = exp(-2*(x - 1).^2);
elseif (ictype ==4) %sinus wave
    ui=sin(2*pi*x);
elseif (ictype==5) %Piecewise continuous
    for i=1:length(x)
        if (x(i) < 0)
            ui(i)=1;
        elseif (0<= x(i) && x(i) < 1)
            ui(i)=1-x(i);
        elseif (x(i) >= 1)
            ui(i)=0;
        end
    end
elseif (ictype==6) %%dirac delta
    for i=1:length(x)
            if ((xend/2) == x(i) )
                ui(i)=10;
            else
                ui(i)=0;
            end
        end
    end
end