function [rhsu] = HeatCRHS1Dsdd(u,time,roota)

% function [rhsu] = HeatCRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D heat equation using central flux

Globals1D;

% Define field differences at faces
drootau = zeros(Nfp*Nfaces,K); drootau(:) = (roota(vmapM).*u(vmapM) - roota(vmapP).*u(vmapP))/2.0;

% impose boundary condition -- Dirichlet BC's
uin  = -u(vmapI); drootau(mapI) = (roota(vmapI).*u(vmapI)-roota(vmapI).*uin)/2.0; 
uout = -u(vmapO); drootau(mapO) = (roota(vmapO).*u(vmapO) - roota(vmapO).*uout)/2.0;

% Compute q and form differences at faces
q = roota.*(rx.*(Dr*u)) - LIFT*(Fscale.*(nx.*drootau));
drootaq = zeros(Nfp*Nfaces,K); drootaq(:)  = (roota(vmapM).*q(vmapM)-roota(vmapP).*q(vmapP))/2.0;

% impose boundary condition -- Neumann BC's
qin  = q(vmapI); drootaq(mapI) = (roota(vmapI).*q(vmapI)- roota(vmapI).*qin )/2.0; 
qout = q(vmapO); drootaq(mapO) = (roota(vmapO).*q(vmapO)- roota(vmapO).*qout)/2.0;

% compute right hand sides of the semi-discrete PDE
rhsu = rx.*(Dr*(roota.*q)) - LIFT*(Fscale.*(nx.*drootaq));
return
