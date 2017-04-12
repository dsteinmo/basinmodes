function [IM] = InterpMatrix1D(rout)

% function [IM] = InterpMatrix1D
% Purpose: Compute local elemental interpolation matrix
% that maps data on the current element to the nodes in rout

Globals1D;
 
% compute Vandermonde at (rout,sout)
Vout = Vandermonde1D(N, rout);

% build interpolation matrix
IM = Vout*invV;
return

  
  
