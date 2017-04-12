%works, but do algebra to reduce number of computations
%Also, problems in lower left corner at low N. 
%find out why, check your 1D version where you had problems
%and fixed it.

%doing this seemed to make the problem in the bottom left corner worse.
%did I make a mistake, or did I correct stuff that brought the problem
%to light more?

%on bc's: see if you can impose them just by modifying
%+ traces.
function [rhsu,Mu] = HelmRHS2DsddLDGv2(u,roota)

% function [rhsu] = PoissonRHS2D(u)
% Purpose  : Evaluate RHS flux in 2D Heat/Poisson equation using a stabilized internal penalty flux

Globals2D;

updownflag=1; %must be either +/-1H

% Define field differences at faces and impose Dirichlet BCs
% Note that Neumann BC's require no work (Natural)
%du = zeros(Nfp*Nfaces,K);
drootau = zeros(Nfp*Nfaces,K); 


rootauM = zeros(Nfp*Nfaces,K);
rootauP = zeros(Nfp*Nfaces,K);

rootauM(:) = roota(vmapM).*u(vmapM);
rootauP(:) = roota(vmapP).*u(vmapP);
%modify + traces to impose BC's (only in dirichet case)
rootauP(mapD) = -rootauM(mapD);

%define field differences
drootau(:) = rootauM - rootauP;
%drootau(mapD) = 2*roota(vmapD).*u(vmapD); what we used to do for BC's

% Compute qx and qy, define differences and impose Neumann BC's
[dudx,dudy] = Grad2D(u);

% Compute DG gradient with central fluxes
% i.e. compute x & y components of q
%fluxxu = nx.*du/2.0; qx = dudx - LIFT*(Fscale.*fluxxu); 
%fluxyu = ny.*du/2.0; qy = dudy - LIFT*(Fscale.*fluxyu);


rootauavg = (rootauM + rootauP)/2;
rootaustarx = rootauavg + updownflag*nx.*drootau;
rootaustary = rootauavg + updownflag*ny.*drootau;

rootaustarx(mapD) = rootauavg(mapD) - updownflag*nx(mapD).*drootau(mapD);
rootaustary(mapD) = rootauavg(mapD) - updownflag*ny(mapD).*drootau(mapD);

fluxxu = nx.*(rootauM - rootaustarx); qx = roota.*dudx - LIFT*(Fscale.*fluxxu); 
fluxyu = ny.*(rootauM - rootaustary); qy = roota.*dudy - LIFT*(Fscale.*fluxyu);

% Compute minimum height of elements either side of each edge
hmin = min(2*J(vmapP)./sJ(mapP), 2*J(vmapM)./sJ(mapM));
tau = reshape(Np./hmin, Nfp*Nfaces, K); 

%now the u stuff

rootaqxM = zeros(Nfp*Nfaces,K); rootaqxP = zeros(Nfp*Nfaces,K);
rootaqyM = zeros(Nfp*Nfaces,K); rootaqyP = zeros(Nfp*Nfaces,K);

rootaqxM(:) = roota(vmapM).*qx(vmapM); rootaqxP(:) = roota(vmapP).*qx(vmapP);
rootaqyM(:) = roota(vmapM).*qy(vmapM); rootaqyP(:) = roota(vmapP).*qy(vmapP);
%modify + traces to impose BC's (only in Neuman case)
rootaqxP(mapN) = -rootaqxM(mapN);
rootaqyP(mapN) = -rootaqyM(mapN);

% Evaluate jumps in components of q at element interfaces
drootaqx = zeros(Nfp*Nfaces,K);
%drootaqx(:) = roota(vmapM).*qx(vmapM) - roota(vmapP).*qx(vmapP);
drootaqx = rootaqxM - rootaqxP;
%drootaqx(mapN) = 2*roota(vmapN).*qx(vmapN); %Need to do this in case BC's on u are Neumann

drootaqy = zeros(Nfp*Nfaces,K);
drootaqy = rootaqyM - rootaqyP;
%drootaqy(:) = roota(vmapM).*qy(vmapM) - roota(vmapP).*qy(vmapP);
%drootaqy(mapN) = 2*roota(vmapN).*qy(vmapN); %Need to do this in case BC's on u are Neumann

% Evaluate flux function

rootaqxavg = (rootaqxM + rootaqxP)/2;
rootaqyavg = (rootaqyM + rootaqyP)/2;
rootaqstarx = rootaqxavg - updownflag*nx.*drootaqx;
rootaqstary = rootaqyavg - updownflag*ny.*drootaqy;

rootaqstarx(mapN) = rootaqxavg(mapN) + updownflag*nx(mapN).*drootaqx(mapN);
rootaqstary(mapN) = rootaqxavg(mapN) + updownflag*ny(mapN).*drootaqy(mapN);


fluxq = (nx.*(rootaqxM-rootaqstarx) + ny.*(rootaqyM-rootaqstary) + tau.*drootau/2);

% Compute right hand side
%divq = Div2D(qx,qy);
divrootaq = Div2D(roota.*qx,roota.*qy);

% compute right hand side residual
%rhsu = J.*((invV'*invV)*(divq - LIFT*(Fscale.*fluxq)));
rhsu = J.*((invV'*invV)*(divrootaq - u - LIFT*(Fscale.*fluxq)));
Mu = J.*((invV'*invV)*u);
return;
