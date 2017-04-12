% Driver script for solving the 2D Poisson equation
% Note that if you use Neumann conditions the solution will have some
% unknown constant added to it (non-uniqueness of potential function)
% Spatially dependent diffusivity working. 

Globals2D;

% Polynomial order used for approximation 
N = 4;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('circA01.neu'); %Dirichlet circle (no curvilinear elements)
%[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('circA01derek.neu'); %same, but with neumann conditions instead of dirichlet.
%[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('circA01a.neu'); %diamond thing.
%[Nv,VX,VY,K,EToV] = MeshReaderGambitBC2D('squareireg2.neu'); %11-element square. [-1,1]x[-1,1]. Dirichlet conditions
%[Nv,VX,VY,K,EToV] = MeshReaderGambitBC2D('squareireg2neu.neu'); %same, but with Neuman conditions

% Initialize solver and construct grid and metric
StartUp2D;

% set up boundary conditions
BuildBCMaps2D;

%set spatially dependent diffusivity
%a = 0*cos(pi*x).*cos(pi*y) + 5;
%a = cos(pi*x).*cos(pi*y);
%a = -2*(pi^2)*sin(pi*x).*sin(pi*y);
a=ones(size(x));
%a = x.*y;

% set up right hand side for homogeneous Poisson 
[A,M] = Poisson2D(); %Setup using PoissonRHS2D.m, doesn't work with inhomogenous BC's. but can do homo neumann or dirichlet
%[A,M] = Poisson2Dsdd(a); %same as above but with spatial dep. diffusivity.
%[A,M] = PoissonIPDG2D_derek(); % Setup using PoissonIPDG2D.m  %this is LDG?, works with inhomogenous BC's

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
uD(mapD) = sin(pi*Fx(mapD)).*sin(pi*Fy(mapD));



% set up Neumann boundary conditions (this only kicks in if mesh supports it)
qN = zeros(Nfp*Nfaces, K);
qN(mapN) = nx(mapN).*(pi*cos(pi*Fx(mapN)).*sin(pi*Fy(mapN))) + ...
           ny(mapN).*(pi*sin(pi*Fx(mapN)).*cos(pi*Fy(mapN))) ;


% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D(uD, qN);


% set up right hand side forcing
force = -2*(pi^2)*cos(pi*x).*cos(pi*y);  %use this for a = 1
%force = -2*(pi^2)*sin(pi*x).*sin(pi*y);
%force = pi*y.*cos(pi*x).*sin(pi*y) - 2*pi^2*x.*y.*sin(pi*x).*sin(pi*y) + pi*x.*sin(pi*x).*cos(pi*y);  %use this for a=xy

%rhs = -MassMatrix*(J.*force) + Aqbc;  %Use this if you're solvling -Lap*u= f
rhs = MassMatrix*(J.*force) - Aqbc;  %Use this if you're solving Lap*u = f

% solve system
%Matlab uses UMFPACK's symbolic LU with automatic re-ordering, then numeric
%LU, followed by UMFPACK's triangular solve
u = A\rhs(:);

u = reshape(u, Np, K);

uexact = sin(pi*x).*sin(pi*y);  %exact solution
%uexact = cos(pi*x).*cos(pi*y);

%u = u-1.08342535e9 -7.8 ;
%u = u - 1.5234e9;
%u = u -1.767039695e11 - 3.9697e6;

 figure(1);
 PlotField2D(N,x,y,u); view([0 90]); colorbar;
 title('Solution - homogenous Neumann and Central flux');
 title('DG Solution');
 
 return;
% figure(2);
% PlotMesh2D;
% figure(3);
% PlotField2D(N,x,y,force); view([0 90]); colorbar;
% title('forcing');
% figure(4);
% PlotField2D(N,x,y,a); view([0 90]); colorbar;
% title('diffusivity');
% figure(5);
% PlotField2D(N,x,y,uexact); view([0 90]); colorbar;
% title('Exact Solution');
% error = norm(uexact(:)-u(:),2)/norm(uexact(:),2)


%look at eigenvalues
% numeigs=21;
% 
% [V,d]=eigs(A,numeigs,'SM');   
% %if using 'eig'
% % [V,d] = eig(full(A));
% %    d=diag(d);
% %    [d,inds] = sort(d,'descend');
% %    V = V(:,inds);
% %end if using 'eig
% phi = cell(numeigs,1);
% for jj=1:numeigs
%     phi{jj} = V(:,jj);
%     phi{jj} = reshape(phi{jj},Np,K);
%     figure(jj);
%     PlotField2D(N,x,y,phi{jj}); view([0 90]); colorbar;
% end
