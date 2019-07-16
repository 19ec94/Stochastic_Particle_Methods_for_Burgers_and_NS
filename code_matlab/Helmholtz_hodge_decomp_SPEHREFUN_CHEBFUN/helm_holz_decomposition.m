%%
%f = spherefunv( @(x,y,z) y.*z.*cos(x.*y.*z), ...
%                @(x,y,z) x.*z.*sin(4*x+.1*y+5*z.^2), @(x,y,z) 1+x.*y.*z );
f = spherefunv( @(x,y,z) y.*z.*cos(x.*y.*z), ...
                @(x,y,z) x.*z.*sin(4*x+.1*y+5*z.^2), @(x,y,z) 0);
quiver3( f ), view([-36,8])

%%
f = tangent( f );
quiver3( f ), view([-36 8])
%%
phi = spherefun.poisson(divergence(f),0,251);
quiver3(gradient(phi)),hold on,
contour(phi,'b-'), view([-36,8]),hold off
%%
psi = spherefun.poisson(vorticity(f),0,251);
quiver3(curl(psi)),hold on,
contour(psi,'r-')
view([-36,8]),hold off
%%
subplot(1,3,1)
quiver3( gradient( phi ) ), title('Curl-free'), view([-36 8])
subplot(1,3,2)
quiver3( curl( psi ) ), title('Divergence-free'), view([-36 8])
subplot(1,3,3)
quiver3( f ), title('Tangent vector field'), view([-36 8])
%%