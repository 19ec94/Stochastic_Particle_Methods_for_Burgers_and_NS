xend = 2;
tend = 2;
N=10;
dx = xend/N;
dt = 0.01;
x = 0:dx:xend;
nt = floor(tend/dt);
dt = tend / nt;
ictype=3
u0 = uinit(x,ictype);
u = u0;
unew = 0*u;
for i = 1 : nt,
        D=0.01;
        fminus = 0.5*( f(u(2:end-1)) + f(u(1:end-2)) );
        fplus = 0.5*( f(u(2:end-1)) + f(u(3:end)) );
        unew(2:end-1) = u(2:end-1) + dt*(D*(u(3:end)- ...
            2*u(2:end-1)+u(1:end-2))/(dx)^2 ...
            - (fplus-fminus)/dx );
        unew(1) = u(1);
        unew(end) = u(end);
end
u = unew;
U(i,:) = u(:);

U=[u0;U];
T=0:dt:tend;

figure(1)
surf(x,T,U)
shading interp
xlabel('x'), ylabel('t'), zlabel ('u(x,t)');
grid on
colormap('Gray');
