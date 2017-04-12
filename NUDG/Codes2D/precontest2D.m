
% Driver script for solving the 2D Helmholtz equation with spatially
% varying diffusivity

Globals2D;

%read in mesh.
load('boxnocornersK288.mat');


%build coarse grid and operator, store the factors.
Ncoarse=4;
Nfine=6;
N=Ncoarse;

%Initialize coarse grid and required maps
StartUp2D; BuildBCMaps2D;
Npcoarse = Np; %store number of points per element on coarse grid
xcoarse=x; ycoarse=y;
rcoarse=r; scoarse=s; %save coarse distribution of Warbuton's poitns for later
acoarse=x.*y+2;
disp('building coarse operator...');
%tic; Acoarse = Helm2DsddLDG(acoarse); toc;
tic; Acoarse = Helm2Dsdd(acoarse); toc;
disp('done.');
disp('factorizing coarse operator...');
tic; [ll,uu,pp,qq] = lu(Acoarse); toc;
disp('done');
%done factorizing coarse operator

%Now get the fine distribution of Warbuton's points
[xfine,yfine] = Nodes2D(Nfine); 
[rfine,sfine] = xytors(xfine,yfine);

%Construct interpolation matrix (coarse -> fine)
InterpMat = InterpMatrix2D(rfine,sfine);

%Now build real grid and operators

% Polynomial order used for approximation 
N = Nfine;
StartUp2D; BuildBCMaps2D;

%Construct coarsening matrix (fine -> coarse)
CoarseMat = InterpMatrix2D(rcoarse,scoarse);


%convergence problems near corners on a box-domain
%probably because the normal changes discontinuously
%within a single element

%oddly, the divergence away from the correct solution
%always seems to happen at a lone corner with both flux choices.


%note that central flux is going to have bad convergence
%for odd choices of N.
 
      

%possible remedies: (1) add more elements near the corner
%                   (2) increase order so stencil talks to more elements
%                   than just nearby ones.

%Things to check: 
%There may be a bug in H&W's implementation of central flux
%so write it "the long way" and compare

%still could do algebraic simplification to your LDG implementation
%call it, "ldgOPT" or something.




%set spatially dependent diffusivity
%a = 0*cos(pi*x).*cos(pi*y) + 5;
%a = cos(pi*x).*cos(pi*y);
%a = -2*(pi^2)*sin(pi*x).*sin(pi*y);
%a=ones(size(x));
a = x.*y+2;
roota = sqrt(a); %for the operator.

%a = x.*y; %Note this may be a bad choice: singular at x,y=0,0

% set up right hand side for homogeneous Helmholtz problem 
% disp('building operator...');
% tic;
% [A,M] = Helm2Dsdd(a); %same as above but with spatial dep. diffusivity.
% toc;
% disp('done.');

% set up right hand side for homogeneous Helmholtz problem with LDG flux
%[ALDG,MLDG] = Helm2DsddLDG(a); %same as above but with spatial dep. diffusivity.

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
%uD(mapD) = sin(pi*Fx(mapD)).*sin(pi*Fy(mapD));

% set up Neumann boundary conditions (this only kicks in if mesh supports it)
qN = zeros(Nfp*Nfaces, K);
%qN(mapN) = nx(mapN).*(pi*cos(pi*Fx(mapN)).*sin(pi*Fy(mapN))) + ...
%          ny(mapN).*(pi*sin(pi*Fx(mapN)).*cos(pi*Fy(mapN))) ;


% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D(uD, qN);

%exact solution for Neumann conditions
uexact = cos(pi*x).*cos(pi*y);  

% set up right hand side forcing
%force = pi*y.*cos(pi*x).*sin(pi*y) - 2*pi^2*x.*y.*sin(pi*x).*sin(pi*y) + pi*x.*sin(pi*x).*cos(pi*y) - sin(pi*x).*sin(pi*y);  %use this for a=xy, Dirichlet
force1 = -pi*y.*sin(pi*x).*cos(pi*y)-2*pi^2*x.*y.*cos(pi*x).*cos(pi*y)-pi*x.*cos(pi*x).*sin(pi*y)-cos(pi*x).*cos(pi*y); %use this for a=xy, Neumann

%what if force = Divergence of \vec{stuff}?
%set components of the 2-vector
stuffx = -pi*(x.*y+2).*sin(pi*x).*cos(pi*y)-(1/2/pi)*sin(pi*x).*cos(pi*y);
stuffy = -pi*(x.*y+2).*cos(pi*x).*sin(pi*y)-(1/2/pi)*cos(pi*x).*sin(pi*y);
divstuff = Div2D(stuffx,stuffy);  %compute weak divergence
dstuffx = zeros(Nfp*Nfaces,K); dstuffy = zeros(Nfp*Nfaces,K);
dstuffx(:) = (stuffx(vmapM)-stuffx(vmapP));  %form field differences
dstuffy(:) = (stuffy(vmapM)-stuffy(vmapP));
fluxstuff = (nx.*dstuffx + ny.*dstuffy)/2.0; %compute central flux
force = divstuff - LIFT*(Fscale.*fluxstuff);
%force = J.*((invV'*invV)*(divstuff - LIFT*(Fscale.*fluxstuff)));

%figure(18);
%PlotField2D(N,x,y,force); shading interp; colorbar; view([0 90]);
%title('RHS'); drawnow;

%forcecoarse = CoarseMat*force;

%doesn't work because of the globals nonsense.
% figure(19);
% PlotField2D(Ncoarse,xcoarse,ycoarse,forcecoarse); shading interp; colorbar; view([0 90]);
% title('RHS - coarsened'); drawnow;


%rhs = -MassMatrix*(J.*force) + Aqbc;  %Use this if you're solvling -Lap*u= f
rhs = MassMatrix*(J.*force) + Aqbc;  %Use this if you're solving Lap*u = f
%rhs = MassMatrix*(J.*force);

% solve system
%Matlab uses UMFPACK's symbolic LU with automatic re-ordering, then numeric
%LU, followed by UMFPACK's triangular solve
% disp('inverting with backslash...');
% tic;
% u = A\rhs(:);
% toc;
% u = reshape(u, Np, K);
% 
% 
% disp('inverting with gmres (no preconditioner)...');
% tic;
% ugmres = gmres(@(foo)HelmRHS2DsddLDGvec(foo,roota),rhs(:),[],1e-6,Np*K);
% toc;
% ugmres = reshape(ugmres,Np,K);
% 
% disp('inverting with gmres (preconditioner)...');
% tic;
% ugmres2 = gmres(@(foo)HelmRHS2Dsddvec(foo,roota),rhs(:),[],1e-6,Np*K,@(bar)mypreconDG2D(bar,ll,uu,pp,qq,InterpMat,CoarseMat,Np,Npcoarse,K));
% toc;
% ugmres2 = reshape(ugmres2,Np,K);

disp('inverting with gmres (preconditioner) and central fine flux...');
tic;
ugmres2 = gmres(@(foo)HelmRHS2Dsddvec(foo,roota),rhs(:),[],1e-6,Np*K,@(bar)mypreconDG2D(bar,ll,uu,pp,qq,InterpMat,CoarseMat,Np,Npcoarse,K));
toc;
ugmres2 = reshape(ugmres2,Np,K);

disp('inverting with 2-lvl solver...');
tic;
mytol = 1e-3;
cvgce = inf;
uguess = ugmres2 + cos(3*pi*x);
%dtau=2.8e-2; %can probably get away with dt=1 on lake-scale
dtau=1e-2;
NumRelax = 100; %number of pseudotimestepping relaxations
It=0;

%main MG loop
while cvgce > mytol
    %relax on u.
    usmoothed = pseudorelax(uguess,dtau,roota,rhs,NumRelax);
  %  figure(4); PlotField2D(N,x,y,usmoothed); view([0 90]); colorbar; title('u_{smoothed}'); drawnow;
    
    
    %compute residual
    ufine = usmoothed;
    resfine = HelmRHS2Dsdd(ufine,roota) - rhs;
    %keyboard;
    %coarsen residual
    rescoarse = CoarseMat*resfine;
    %Now Direct-Solve A*errcoarse = rescoarse
    errcoarse = qq*(uu\(ll\(pp*rescoarse(:))));
    errcoarse = reshape(errcoarse,Npcoarse,K);
    %Interpolate error to the fine grid.
    errfine = InterpMat*errcoarse;
    %correct solution on fine grid
    ufine = ufine + errfine;
    %relax on coarse-grid corrected u.
    ufinesmoothed = pseudorelax(ufine,dtau,roota,rhs,NumRelax);
        
    %calculate cauchy cvgce. criteria
    cvgce = norm(ufinesmoothed-uguess,2)
    if cvgce > 1e8
        disp('main MG loop failed to converge.');
        return;
    end
    
    %rotate for next loop
    uguess = ufinesmoothed;
    It=It+1
end

toc;
figure(4); PlotField2D(N,x,y,ufinesmoothed); view([0 90]); colorbar; title('u_{2lvl}'); drawnow;








figure(1);
PlotField2D(N,x,y,u); view([0 90]); colorbar;
title('u_{central} - backslash'); drawnow;
figure(2);
PlotField2D(N,x,y,ugmres); view([0 90]); colorbar;
title('u_{central} - gmres'); drawnow;
figure(3);
PlotField2D(N,x,y,ugmres2); view([0 90]); colorbar;
title('u_{central} - gmres (preconditioned)'); drawnow;
figure(4); PlotField2D(N,x,y,ufine); view([0 90]); colorbar; title('u_{2lvl}'); drawnow;


relerror = norm(u-uexact,2)/norm(uexact,2)
relerrorgmres = norm(ugmres-uexact,2)/norm(uexact,2)
rellerrorgmres2 = norm(ugmres2-uexact,2)/norm(uexact,2)

return;

%try LU-factorization
% disp('computing LU factorization of central flux...');
% tic;
% [L,U,P,Q] = lu(A);
% toc;
% disp('done.');
% disp('inverting with lu-factors...');
% tic;
% u1 = Q*(U\(L\(P*rhs(:))));
% toc;

disp('Inverting with LDG (backslash)...');
tic;
uLDG = ALDG\rhs(:);
toc;
%norm(uLDG-u1,2)
uLDG = reshape(uLDG, Np, K);


% disp('computing LU factorization of LDG operator...');
% tic;
% [LLDG,ULDG,PLDG,QLDG] = lu(ALDG);
% tic;
% 
% 
% 
% disp('computing symmetric cuthill mckee pre-ordering (central)'); %this is working!
% tic;
% P = symrcm(A);
% AP = A(P,P);
% myrhs = rhs(:); 
% myrhs = myrhs(P);
% toc;
% [Lp,Up,Pp,Qp] = lu(AP);
% disp('inverting with pre-ordering and lu...');
% tic;
% u2 = Qp*(Up\(Lp\(Pp*myrhs)));
% toc;
% u2(P) = u2;  %this is working now
% 
% disp('computing symmetric cuthill mckee pre-ordering (LDG)'); %this is working!
% tic;
% PL = symrcm(ALDG);
% ALDGP = ALDG(PL,PL);
% myrhs = rhs(:); 
% myrhs = myrhs(PL);
% toc;
% [LpLDG,UpLDG,PpLDG,QpLDG] = lu(ALDGP);
% disp('inverting with pre-ordering and lu...');
% tic;
% u2 = QpLDG*(UpLDG\(LpLDG\(PpLDG*myrhs)));
% toc;
% u2(P) = u2;  %this is working now



%uexact = sin(pi*x).*sin(pi*y);  %exact solution for Dirichlet conditions
uexact = cos(pi*x).*cos(pi*y);  %exact solution for Neumann conditions


figure(1);
PlotField2D(N,x,y,u); view([0 90]); colorbar;
%title('Solution - homogenous Neumann and Central flux');
title('u_{central}');
%figure(2);
%PlotMesh2D;
% figure(3);
% PlotField2D(N,x,y,force); view([0 90]); colorbar;
% title('forcing - as a DG divergence');
%figure(4);
%PlotField2D(N,x,y,a); view([0 90]); colorbar;
%title('diffusivity');
%figure(5);
%PlotField2D(N,x,y,uexact); view([0 90]); colorbar; title('exact solution');
centralerror = norm(uexact(:)-u(:),2)/norm(uexact(:),2)
figure(6);
PlotField2D(N,x,y,uLDG); view([0 90]); colorbar;
title('u_{LDG}');
LDGerror = norm(uLDG(:)-uexact(:),2)/norm(uexact(:),2)
% figure(6);
% PlotField2D(N,x,y,force1); view([0 90]); colorbar;
% title('forcing -exact');
% figure(8);
% subplot(4,3,1);
% spy(A); title('Central flux');
% disp(['A Fullness %: ' num2str(nnz(A)/numel(A))]);
% subplot(4,3,2);
% spy(L); title('central L');
% disp(['L Fullness %: ' num2str(nnz(L)/numel(L))]);
% subplot(4,3,3);
% spy(U); title('central U');
% disp(['U Fullness %: ' num2str(nnz(U)/numel(U))]);
% subplot(4,3,4);
% spy(Lp); title('central L (pre-ordered)');
% disp(['L (pre-ordered) Fullness %: ' num2str(nnz(Lp)/numel(Lp))]);
% subplot(4,3,5);
% spy(Up); title('central U (pre-ordered)');
% disp(['U (pre-ordered) Fullness %: ' num2str(nnz(Up)/numel(Up))]);
% subplot(4,3,6);
% spy(AP); title('Central flux (pre-ordered)');
% subplot(4,3,7);
% spy(ALDG); title('LDG flux');
% subplot(4,3,8);
% spy(ALDGP); title('LDG (pre-ordered)');
% subplot(4,3,9);
% spy(LLDG); title('LDG L');
% subplot(4,3,10);
% spy(ULDG); title('LDG U');
% subplot(4,3,11);
% spy(LpLDG); title('LDG L (pre-ordered)');
% subplot(4,3,12);
% spy(UpLDG); title('LDG U (pre-ordered)');