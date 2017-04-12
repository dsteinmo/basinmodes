function [A] = HelmCstab1Dsdd(a)

% function [A] = PoissonCstab1D();
% Purpose: Set up symmetric Poisson matrix with estabilized central fluxes

Globals1D;
A = zeros(K*Np); g = zeros(K*Np,1);

roota = sqrt(a);

% Build matrix -- one column at a time
for i=1:K*Np
    g(i) = 1.0;
    gmat = reshape(g,Np,K);
    %Avec = PoissonCstabRHS1D(gmat);  %stabilized C-flux
    %Avec = PoissonIPstabRHS1D(gmat); %IP flux
    Avec = HelmCstabRHS1Dsdd(gmat,[],roota); %stabilized C-flux with spat. dep. diff.
    A(:,i) = reshape(Avec,K*Np,1);
    g(i)=0.0;
end
return
