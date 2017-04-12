function [A,M] = Helm2Dsddnomass2(a)

% function [A,M] = Poisson2D()
% Purpose: Set up matrix for 2D Poisson equation based on stabilized 
%          internal fluxes on symmetric form (homogenous BC's only)
%a = spatially-dependent diffusivity

roota = sqrt(a);

Globals2D;
g = zeros(K*Np,1);
A = spalloc(K*Np, K*Np, 3*Np);  M = spalloc(K*Np, K*Np, 3*Np); 

% Build matrix -- one column at a time
for i=1:K*Np
    g(i) = 1.0;
    gmat = reshape(g,Np,K);
    [Avec,Mvec] = HelmRHS2Dsddnomass2(gmat,roota);
   
    ids = find(Avec); A(ids,i) = Avec(ids);
    ids = find(Mvec); M(ids,i) = Mvec(ids);
    g(i)=0.0;
end
return
