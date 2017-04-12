function [F] = CutOffFilter1D(Nc,frac)

% function [F] = Filter1D(N,Nc,s)
% Purpose : Initialize 1D filter matrix of size N.
%           Order of exponential filter is (even) s with cutoff at Nc;

Globals1D;

filterdiag = ones(N+1,1);


% Initialize filter function
for i=Nc:N
    filterdiag(i+1) = frac;
end

figure(17);
plot(filterdiag,'.'); drawnow;

F = V*diag(filterdiag)*invV;
return;
