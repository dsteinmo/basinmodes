% Driver script for solving the 1D Poisson equation
Globals1D;

% Polynomial order used for approximation 
N =8;

% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(0,pi,10);

% Initialize solver and construct grid and metric
StartUp1D;

% Set RHS (Mass matrix times forcing)
%force = -(sin(2*x) + sin(x));
%f = J.*((invV'*invV)*force);  %to compare with exact solution

%what if RHS = d/dx of stuff ?
stuff = cos(2*x)/2 + cos(x);
dstuff = zeros(Nfp*Nfaces,K); dstuff(:) = (stuff(vmapM)-stuff(vmapP))/2.0;

f = J.*((invV'*invV)*(rx.*(Dr*stuff) - LIFT*(Fscale.*dstuff) )); %Agrees with exact solution!

%set spatially dependent diffusivity
a = cos(x);

% Set up operator
%[A] = PoissonCstab1D();
[A] = HelmCstab1Dsdd(a);

% Solve Problem
solvec = A\f(:);
u = reshape(solvec,Np,K);

uexact = sin(x);

figure(1);
plot(x(:),u(:),x(:),uexact(:));

norm(u(:)-uexact(:),2)
