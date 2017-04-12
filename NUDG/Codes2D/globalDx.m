function [OP,MM] = globalDx()

% Purpose: Set up the discrete Poisson matrix directly
%          using LDG. The operator is set up in the weak form
%          this guy handles inhomogenous BC's.


% OP*u would give MM*Dx*u. Uses central fluxes (no penalty terms, yet)

Globals2D;

% build local face matrices
massEdge = zeros(Np,Np,Nfaces);
Fm = Fmask(:,1); faceR = r(Fm); 
V1D = Vandermonde1D(N, faceR);  massEdge(Fm,Fm,1) = inv(V1D*V1D');
Fm = Fmask(:,2); faceR = r(Fm); 
V1D = Vandermonde1D(N, faceR);  massEdge(Fm,Fm,2) = inv(V1D*V1D');
Fm = Fmask(:,3); faceS = s(Fm); 
V1D = Vandermonde1D(N, faceS);  massEdge(Fm,Fm,3) = inv(V1D*V1D');

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
MM  = zeros(K*Np*Np, 3);  OP = zeros(K*Np*Np*(1+Nfaces), 3);  

% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)'; 
for k1=1:K 
  if(~mod(k1,1000)) k1, end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

  
  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;   Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  MassMatrixinv = inv(MassMatrix);
  
  OP11 = J(1,k1)*(MassMatrix*Dx);
  
  % Build element-to-element parts of operator
  for f1=1:Nfaces
    k2 = EToE(k1,f1); %get number of neighbouring element along f1
    f2 = EToF(k1,f1); %get the face on element k2 that f1 on k1 talks to.
    

    rows2 = ((k2-1)*Np+1:k2*Np)'*ones(1,Np); 
    cols2 = rows2';  %find the entries of the contributing neighbouring element
    
    fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp);
    vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;  %make '-' trace vertex map and Facemask
    vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;  %make '+' trace vertex map and Facemask
    
    id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
    lnx = nx(id);  lny = ny(id); lsJ = sJ(id);     %get the normal and the Jacobian along face (scalars)
    hinv = max(Fscale(id), Fscale(1+(f2-1)*Nfp, k2)); %compute this for penalty scaling

    Dx2 = rx(1,k2)*Dr + sx(1,k2)*Ds;   %build derivative operators on
    Dy2 = ry(1,k2)*Dr + sy(1,k2)*Ds;   %element k2 (NpxNp matrices)
    
    %Dn1 = lnx*Dx  + lny*Dy ;  %normal derivative for element k1 (along this face)
    %Dn2 = lnx*Dx2 + lny*Dy2;  %normal derivative for element k2 (along this face), both (NpxNp

    mmE = lsJ*massEdge(:,:,f1); %get local edge mass matrix along f1
    
    

    %gtau = 100*2*(N+1)*(N+1)*hinv; % set penalty scaling
    
    switch(BCType(k1,f1))
      case {Dirichlet}    %if we're on a Dirichlet bdry, put in appropriate term
	      %OP11 = OP11 + ( mmE*lnx  ); % had this
          OP11 = OP11 - ( mmE*lnx  ); % subtracted ghost term is negative of '-' term, so just get twice '-' term
          %OP11        = OP11 + 0.5*( gtau*mmE*sqrt(2) - mmE*Dn1*sqrt(2) - Dn1'*mmE*sqrt(2) );
      case {Neuman}
        % nada 
      case {Wall}
        % nada
        otherwise   %if we're not on a boundary ... 
	% interior face variational terms
	%OP11        = OP11 + 0.5*( mmE*lnx );  %had this
    OP11        = OP11 - 0.5*( mmE*lnx ); %trying this... (think this is right, given sign of volume term  above)
    
    % contributions from neighbouring element
    
	OP12 = zeros(Np);
   
    %OP12(:,Fm2) = OP12(:,Fm2) - 0.5*mmE(:,Fm1)*lnx; %had this
    OP12(:,Fm2) = OP12(:,Fm2) + 0.5*mmE(:,Fm1)*lnx; %trying this... (think this is right, given sign of volume term above)
    
    OP(entries(:), :) = [rows1(:), cols2(:), OP12(:)];

    
    %OP12(:,Fm2) =  OP12(:,Fm2) 
    
	entries = entries + Np*Np;

    end 
  end      
  OP(entries(:), :)   = [rows1(:), cols1(:), OP11(:)];
  MM(entriesMM(:), :) = [rows1(:), cols1(:), J(1,k1)*MassMatrix(:)];
  entries = entries + Np*Np; entriesMM = entriesMM + Np*Np;
end  

OP   =   OP(1:max(entries)  -Np*Np,:);  OP   = myspconvert(OP, Np*K, Np*K, 1e-15);
MM   =   MM(1:max(entriesMM)-Np*Np,:);  MM   = myspconvert(MM, Np*K, Np*K, 1e-15);
return
