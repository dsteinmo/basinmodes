function [rhsu] = HelmCstabRHS1Dsdd(u,time,roota)

% function [rhsu] = PoissonCstabRHS1D(u,time)
% Purpose  : Evaluate RHS in 1D Poisson equation on
%            symmetric form using stabilized central flux

Globals1D;

% Define field differences at faces
drootau = zeros(Nfp*Nfaces,K); drootau(:) = (roota(vmapM).*u(vmapM)-roota(vmapP).*u(vmapP))/2.0;

% impose boundary condition -- Dirichlet BC's
%uin  = -u(vmapI); drootau (mapI) = (roota(vmapI).*u(vmapI) - roota(vmapI).*uin )/2.0;
%uout = -u(vmapO); drootau (mapO) = (roota(vmapO).*u(vmapO) - roota(vmapO).*uout)/2.0;

uin  = u(vmapI); drootau (mapI) = (roota(vmapI).*u(vmapI) - roota(vmapI).*uin )/2.0;
uout = u(vmapO); drootau (mapO) = (roota(vmapO).*u(vmapO) - roota(vmapO).*uout)/2.0;

% Compute q
q = roota.*rx.*(Dr*u) - LIFT*(Fscale.*(nx.*drootau));
drootaq = zeros(Nfp*Nfaces,K); drootaq(:) = roota(vmapM).*q(vmapM)-roota(vmapP).*q(vmapP);

% impose boundary condition -- Neumann BC's
% qin  = q(vmapI); drootaq (mapI) = roota(vmapI).*q(vmapI) - roota(vmapI).*qin; 
% qout = q(vmapO); drootaq (mapO) = roota(vmapO).*q(vmapO) - roota(vmapO).*qout;

qin  = -q(vmapI); drootaq (mapI) = roota(vmapI).*q(vmapI) - roota(vmapI).*qin; 
qout = -q(vmapO); drootaq (mapO) = roota(vmapO).*q(vmapO) - roota(vmapO).*qout;

% penalize q-fluxes
tau = 1.0;
fluxq = nx.*(drootaq/2.0 + tau*nx.*drootau);

% compute right hand sides of the semi-discrete PDE
rhsu = J.*((invV'*invV)*(rx.*(Dr*(roota.*q)) - u  - LIFT*(Fscale.*fluxq)));
return
