function [rhsu] = AdvecRHS1D(u,time, a)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
%alpha=1; %central flux
alpha=0; %upwind flux
du = zeros(Nfp*Nfaces,K); 
du(:) = (u(vmapM)-u(vmapP)).*(a*nx(:)-(1-alpha)*abs(a*nx(:)))/2;


% impose boundary condition at x=0
%uin = -sin(a*time);  %inflow condition
%uin = 0;
%du (mapI) = (u(vmapI)- uin ).*(a*nx(mapI)-(1-alpha)*abs(a*nx(mapI)))/2;
%du (mapO) = 0;

% impose periodic BC's (Derek)
du(mapI) = (u(vmapI)-u(vmapO)).*(a*nx(mapI)-(1-alpha)*abs(a*nx(mapI)))/2;
du(mapO) = (u(vmapO)-u(vmapI)).*(a*nx(mapO)-(1-alpha)*abs(a*nx(mapO)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -a*rx.*(Dr*u) + LIFT*(Fscale.*(du));

%try frank's "DSS" (direct stiffness summation) trick
%may need to also lump mass for this to work, or be "rigorous"
%rhsu(vmapM) = 0.5*(rhsu(vmapM)+rhsu(vmapP));
%rhsu(vmapP) = rhsu(vmapM);
return
