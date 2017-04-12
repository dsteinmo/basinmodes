% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation 
N = 4;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,2*pi,20);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
u = sin(x);
%u = sech(2*(x-5));
%u=zeros(Np,K);

% Solve Problem
FinalTime = 10;
[u] = Advec1D(u,FinalTime);
