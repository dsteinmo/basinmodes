% Driver script for solving the 1D Poisson equation
Globals1D;

% Polynomial order used for approximation 
N =8;

% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(0,pi/4,10);

% Initialize solver and construct grid and metric
StartUp1D;

% Set RHS (Mass matrix times forcing)
f = J.*((invV'*invV)*sin(x));  %to compare with exact solution

%set spatially dependent diffusivity
a = cos(x);

% Set up operator
%[A] = PoissonCstab1D();
[A] = PoissonCstab1Dsdd(a);

% Solve Problem
solvec = A\f(:);
u = reshape(solvec,Np,K);

uexact = (pi/4)*log(sec(x) + tan(x))/(log(sqrt(2)+1)) - x;

figure(1);
plot(x(:),u(:),x(:),uexact(:));

norm(u(:)-uexact(:),2)
