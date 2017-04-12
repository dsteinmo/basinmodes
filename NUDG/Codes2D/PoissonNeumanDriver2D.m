%solves Lap u = f with arbitrary neuman boundary conditions
%recovers solvability condition by bordering the matrix 
%as in Trottenburg, or 

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

%split triangles into quarters if you like to get finer mesh.
numhrefines=1;
for j=1:numhrefines
    Hrefine2D(ones(K,1)')
    StartUp2D; BuildBCMaps2D;
end

dof = K*Np;

% set up right hand side for homogeneous Poisson 
[A,M] = PoissonIPDG2D(); % Setup using PoissonIPDG2D.m

%choose which method to deal with solvability condition
%method = 'singvectorremove';
method = 'bordering trick';
%method = 'none';


OP = [A ones(Np*K,1);
      ones(1,Np*K) 0]; %lower border should actually have
                       %quadrature weights, but this seems to work.

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
uD(mapD) = sin(pi*Fx(mapD)).*sin(pi*Fy(mapD));

% set up Neumann boundary conditions
% dg uses auxiliary variable 'q = grad u'
% so dq/dn = grad u dot n
qN = zeros(Nfp*Nfaces, K);
qN(mapN) = nx(mapN).*(pi*cos(pi*Fx(mapN)).*sin(pi*Fy(mapN))) + ...
           ny(mapN).*(pi*sin(pi*Fx(mapN)).*cos(pi*Fy(mapN))) + randn(size(mapN))*1e-1;

% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D(uD, qN);


uexact = sin(pi*x).*sin(pi*y);
rhs = -2*pi.^2*sin(pi*x).*sin(pi*y);

% set up right hand side forcing

rhs = -MassMatrix*(J.*rhs) + Aqbc;

%% solve system
if strcmp(method, 'singvectorremove')
    %remove null-singular vector
    [u0,sig] = SVDS(A,1,0);    
    %%this works:
    %mymat = u0*u0';
    %u = A\(rhs(:)-mymat*rhs(:));
    
    
    %shortcut - this works too.
    u0_squared = mean(u0)^2;
    u = A\(rhs(:) - u0_squared*sum(rhs(:)));
    
    
    
elseif strcmp(method,'bordering trick')
    u = OP\[rhs(:); 0];
    u = u(1:end-1);
else
    u = A\rhs(:);
end

u = reshape(u, Np, K);

%close all;
figure(1); clf;
subplot(3,2,1);
PlotField2D(N,x,y,u); shading interp; colorbar; view([0 90]); axis tight; drawnow;
title('DG solution');
subplot(3,2,2);
PlotField2D(N,x,y,uexact); shading interp; colorbar; view([0 90]); axis tight; drawnow;
title('Exact Solution');


%get gradients
[ux,uy] = Grad2D(u);
[uexactx,uexacty] = Grad2D(uexact);

subplot(3,2,3);
PlotField2D(N,x,y,ux); shading interp; colorbar; view([0 90]); axis tight; drawnow;
title('DG u_x');
subplot(3,2,4);
PlotField2D(N,x,y,uy); shading interp; colorbar; view([0 90]); axis tight; drawnow;
title('DG u_y');

subplot(3,2,5);
PlotField2D(N,x,y,uexactx); shading interp; colorbar; view([0 90]); axis tight; drawnow;
title('exact u_x');
subplot(3,2,6);
PlotField2D(N,x,y,uexacty); shading interp; colorbar; view([0 90]); axis tight; drawnow;
title('exact u_y');

figure(2);
PlotMesh2D; axis tight;
title('FEM mesh');

[V,d] = eigs(OP,6,'SM');

d(:)

figure(12);
pf2d(N,x,y,reshape(V(1:end-1,3),Np,K));
