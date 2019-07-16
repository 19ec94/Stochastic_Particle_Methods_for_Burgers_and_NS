%% GRID setup 
NWire=1;
I=10;
xW=0.00;
yW=0.00;
minZ=-5; maxZ=5;
nZ=4101; %I am going to discretize the wite into 4100 pieces 
radiusW=0.001;

N=401; %400 cells
L=0.005*2; %length of the wire is 0.01
minX= -L; maxX=L; %minX=-0.01, maxX=0.01;
minY= -L; maxY=L; %minX=-0.01, maxX=0.01;
minR=0.001;
%% Constants
mu0=4*pi*1e-7;
K= I*(mu0/(4*pi));
zW = linspace(minZ,maxZ,nZ);
dLz=zW(2)-zW(1); 
dLx=0; dLy=0;
dl = [dLx dLy dLz]; %only z component, similar to vorticity
x= linspace(minX,maxX,N);
y= linspace(minY,maxY,N);
dx = x(2)-x(1); dy=y(2)-y(1);

[xG, yG] = meshgrid(x,y); %create a grid in x-y plane 
RG = sqrt(xG.^2+yG.^2); % length of a distance vector from origin
%%
Bx = zeros(N,N);
By = zeros(N,N);
Bz = zeros(N,N);

for m=1:nZ
    Rx = xG -xW;
    Ry = yG -yW;
    Rz = 0  -zW(m);

    R = sqrt(Rx.^2+Ry.^2+Rz.^2);
    R3 =R.^3;
    
    Bx = Bx + K .* (dLy .* Rz- dLz .* Ry) ./R3;
    By = By + K .* (dLz .* Rz- dLx .* Rz) ./R3;
    Bz = Bz + K .* (dLx .* Rz- dLy .* Rx) ./R3;
end
B = sqrt(Bx.^2 +By.^2 +Bz.^2 );
%%
surf(xG,yG,B)




