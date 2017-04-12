clear; close all;
% Driver script for solving the 2D Poisson equation
Globals2D;

% Polynomial order used for approximation 
N = 4;

% Read in Mesh
%[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('circA01.neu'); %dirichlet
[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('circA01derek.neu'); %neuman




% Initialize solver and construct grid and metric
StartUp2D;

% set up boundary conditions
BuildBCMaps2D;

%hrefine a few times if you like.
numhrefines=0;
for j=1:numhrefines
    Hrefine2D(ones(K,1)')
    StartUp2D; BuildBCMaps2D;
end



dof = K*(Np);

%a=ones(size(x));
a=x.*y+2;

% set up right hand side for homogeneous Poisson 
%[A,M] = Poisson2D(); % Setup using PoissonRHS2D.m
[A,M] = PoissonIPDG2D(); % Setup using PoissonIPDG2D.m

%do bordering trick on the operator for poisson-neuman
BORDERING_TRICK =true;
OP = [A ones(Np*K,1);
      ones(1,Np*K) 0];

tic;
%[A,M] = PoissonIPDG2Dsdd(a);
toc;

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
uD(mapD) = sin(pi*Fx(mapD)).*sin(pi*Fy(mapD));

% set up Neumann boundary conditions
qN = zeros(Nfp*Nfaces, K);
qN(mapN) = nx(mapN).*(pi*cos(pi*Fx(mapN)).*sin(pi*Fy(mapN))) + ...
           ny(mapN).*(pi*sin(pi*Fx(mapN)).*cos(pi*Fy(mapN))) ;

% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D(uD, qN);

[thetam,rm] = cart2pol(x,y);

uexact = sin(pi*x).*sin(pi*y);
%uexact = cos(pi*rm/2);
%rhs = (1./rm).*(-pi/2*sin(pi*rm/2)-pi^2/4*rm.*cos(pi*rm/2));
%rhs= (1./rm).*(-0.5*(rm.^2.*cos(thetam).*sin(thetam)+2).*sin(pi*rm/2) - rm.^2.*cos(thetam).*sin(thetam).*sin(pi*rm/2)*pi - 0.25*rm.*(rm.^2.*cos(thetam).*sin(thetam)+2).*cos(pi*rm/2)*pi^2);
rhs = -2*pi.^2*sin(pi*x).*sin(pi*y);

% set up right hand side forcing

%rhs = -2*(pi^2)*sin(pi*x).*sin(pi*y);
rhs = -MassMatrix*(J.*rhs) + Aqbc;

% solve system



if BORDERING_TRICK == true
    u = OP\[rhs(:); 0];
    u = u(1:end-1);
else
    u = A\rhs(:);
end

u = reshape(u, Np, K);


figure(1);
PlotField2D(N,x,y,u); shading interp; colorbar; view([0 90]); axis tight; drawnow;
title('DG solution');
figure(2);
PlotField2D(N,x,y,uexact); shading interp; colorbar; view([0 90]); axis tight; drawnow;
figure(3);
PlotMesh2D; axis tight;

%agreement won't be amazing,
%since uexact has non-zero mean in this case.

norm(u(:)-uexact(:),2)/norm(u(:),2)