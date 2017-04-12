function [F] = Filter2D(Norder,Nc,sp)

% function [F] = Filter2D(Norder,Nc,sp)
% Purpose : Initialize 2D filter matrix of order sp and cutoff Nc

Globals2D;

filterdiag = ones((Norder+1)*(Norder+2)/2,1);
alpha = -log(eps);

% build exponential filter
sk = 1;
for i=0:Norder
  for j=0:Norder-i
    if (i+j>=Nc)
      filterdiag(sk) = exp(-alpha*((i+j - Nc)/(Norder-Nc))^sp);
    end
    sk = sk+1;
  end
end
% figure(18);
% plot(filterdiag,'.');
% drawnow;
F = V*diag(filterdiag)*invV;
return;
