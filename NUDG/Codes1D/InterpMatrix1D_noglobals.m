function [IM] = InterpMatrix1D_noglobals(rout,N,invV)

% function [IM] = InterpMatrix1D
% Purpose: Compute local elemental interpolation matrix
% that maps data on the current element to the nodes in rout
 
% compute Vandermonde at (rout,sout)
Vout = Vandermonde1D(N, rout);

% build interpolation matrix
IM = Vout*invV;
return

  
  
