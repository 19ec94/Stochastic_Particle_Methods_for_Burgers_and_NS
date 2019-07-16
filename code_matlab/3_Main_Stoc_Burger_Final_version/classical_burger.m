xstart=0;
xend = 1; % x-axis size.
tend = 0.5; % t-axis size.
N = 100;
x = linspace(xstart,xend,N+1);
dx = x(2)-x(1); % Grid spacing
dt = 0.001;
%nt = floor(tend/dt);
%dt = tend / nt;
%Set up the initial solution values.
%u0 = uinit(x,ictype,xend); %Call to the function "uinit".
u0 =zeros(1,N);
for i=1:N
    u0(i)=sin(2*pi* (x(i)+(x(i+1)-x(i))/2));
    %u0(1,i)= exp(-2*((x(1,i)+(x(1,i+1)-x(1,i))/2) - 1).^2);
    %if (i ~= N/2)
    %    u0(1,i)=0;
    %else
    %    u0(1,i)=10;
    %end
end
for i=1:N
    x_new(i)= x(i)+(x(i+1)-x(i))/2;
end
figure(1)
hold on
plot(x_new,u0)
u = u0;
unew = 0*u;
for i = 0:dt:tend
    D=0.01;
    fminus = 0.5*( f(u(2:end-1)) + f(u(1:end-2)) );
    fplus = 0.5*( f(u(2:end-1)) + f(u(3:end)) );
    unew(2:end-1) = u(2:end-1) + dt*(D*(u(3:end)-2*u(2:end-1)+u(1:end-2))/(dx)^2 ...
        - (fplus-fminus)/dx );
    unew(1) = u(1);
    unew(end) = u(end);
    u = unew;
end


figure(1)
hold on
plot(x_new,u)