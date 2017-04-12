function [F] = Filter1D(N,Nc,s,V,invV)

% function [F] = Filter1D(N,Nc,s)
% Purpose : Initialize 1D filter matrix of size N.
%           Order of exponential filter is (even) s with cutoff at Nc;

filterdiag = ones(N+1,1);
alpha = -log(eps);

% Initialize filter function
for i=Nc:N
    filterdiag(i+1) = exp(-alpha*((i-Nc)/(N-Nc))^s);
end;

figure(17);
plot(filterdiag,'.'); drawnow;

F = V*diag(filterdiag)*invV;
return;
