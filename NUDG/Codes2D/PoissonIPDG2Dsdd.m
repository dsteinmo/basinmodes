function [OP,MM] = PoissonIPDG2Dsdd(a)

% Purpose: Set up the discrete Poisson matrix directly
%          using LDG. The operator is set up in the weak form
%          this guy handles inhomogenous BC's.

%Note to self: Will look a little funny because it uses weak form and NOT
%strong form...


%I think this in fact uses an internal penalty flux so it's easier to 
%not deal with 'qx,qy' auxiliary variables in the interior face variational terms.

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

%calculate 'root a' (square root of diffusivity)
rta = sqrt(a);


% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)'; 
for k1=1:K 
  if(~mod(k1,1000)) k1, end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

  %get local 'root a' operator and write as a diagonal matrix
  rtak1 = rta(:,k1);
  rtak1mat = diag(rtak1);
  ak1mat = diag(a(:,k1));
  
  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;   Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

  MassMatrixinv = inv(MassMatrix);
  
  OP11 = J(1,k1)*(Dx'*MassMatrix*rtak1mat*MassMatrixinv*rtak1mat*MassMatrix*Dx + ...
                  Dy'*MassMatrix*rtak1mat*MassMatrixinv*rtak1mat*MassMatrix*Dy);
  
  %OP11 = 2*J(1,k1)*(Dx'*MassMatrix*Dx + Dy'*MassMatrix*Dy);  %volume contribution to the laplacian.
                                                           %...transpose?
                                                           %Yes, because operator is stored in
                                                           %column vector temporarily.
  % Build element-to-element parts of operator
  for f1=1:Nfaces
    k2 = EToE(k1,f1); %get number of neighbouring element along f1
    f2 = EToF(k1,f1); %get the face on element k2 that f1 on k1 talks to.
    
    %get k2's 'root a' operator
    rtak2 = rta(:,k2);
    rtak2mat = diag(rtak2);

    rows2 = ((k2-1)*Np+1:k2*Np)'*ones(1,Np); 
    cols2 = rows2';  %find the entries we're contributing to (for k2) in the global operator
    
    fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp);
    vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;  %make '-' trace vertex map and Facemask
    vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;  %make '+' trace vertex map and Facemask
    
    id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
    lnx = nx(id);  lny = ny(id); lsJ = sJ(id);     %get the normal and the Jacobian along face (scalars)
    hinv = max(Fscale(id), Fscale(1+(f2-1)*Nfp, k2)); %compute this for penalty scaling

    Dx2 = rx(1,k2)*Dr + sx(1,k2)*Ds;   %build derivative operators on
    Dy2 = ry(1,k2)*Dr + sy(1,k2)*Ds;   %element k2 (NpxNp matrices)
    
    Dn1 = lnx*Dx  + lny*Dy ;  %normal derivative for element k1 (along this face)
    Dn2 = lnx*Dx2 + lny*Dy2;  %normal derivative for element k2 (along this face), both (NpxNp)

    mmE = lsJ*massEdge(:,:,f1); %get local edge mass matrix along f1
    
    

    gtau = 100*2*(N+1)*(N+1)*hinv; % set penalty scaling
    
    switch(BCType(k1,f1))
      case {Dirichlet}    %if we're on a Dirichlet bdry, put in appropriate term
	      OP11 = OP11 + ( gtau*mmE*rtak1mat - mmE*rtak1mat*Dn1 - Dn1'*MassMatrix*rtak1mat*MassMatrixinv*rtak1mat*mmE ); % ok
          %OP11        = OP11 + 0.5*( gtau*mmE*sqrt(2) - mmE*Dn1*sqrt(2) - Dn1'*mmE*sqrt(2) );
      case {Neuman}
        % nada 
      case {Wall}
        % nada
        otherwise   %if we're not on a boundary ... 
	% interior face variational terms
	OP11        = OP11 + 0.5*( gtau*mmE*rtak1mat - mmE*rtak1mat*Dn1 - Dn1'*MassMatrix*rtak1mat*MassMatrixinv*rtak1mat*mmE );
    %OP11        = OP11 + 0.5*( gtau*mmE*sqrt(2) - mmE*Dn1*sqrt(2) - Dn1'*mmE*sqrt(2) );
    % contributions from neighbouring element
    
	OP12 = zeros(Np);
    
    
	OP12(:,Fm2) =             - 0.5*( gtau*mmE(:,Fm1)*rtak1mat(Fm1,Fm1) ); %my version
	OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*rtak1mat(Fm1,Fm1)*Dn2(Fm2,:) );
    
    %OP12(:,Fm2) =             - 0.5*( gtau*mmE(:,Fm1)*sqrt(2) );  %a=2 original code
	%OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:)*sqrt(2) ); %a=2original code
    
    %OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1)*2 ); %a=2, original code.
	
    %my attempt
    OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*MassMatrix*rtak1mat*MassMatrixinv*rtak1mat*mmE(:, Fm1) );
    OP(entries(:), :) = [rows1(:), cols2(:), OP12(:)];

    
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
