function [OP,MM] = CurvedPoissonIPDG2Dsdd(a)

% function [OP,MM] = CurvedPoissonIPDG2D()
% Purpose: Set up a discrete Poisson matrix and mass matrix
%          using Interior Penalty Discontinuous Galerkin (IPDG).
 
Globals2D;

NGauss = gauss.NGauss;

% build DG derivative matrices
MM  = zeros(K*Np*Np, 3);  OP = zeros(K*Np*Np*(1+Nfaces), 3); 

%Derek - calculate 'root a' (square root of diffusivity)
roota = sqrt(a);
%Derek - interpolate to cubature nodes
cRoota = cub.V*roota;
cubA = cub.V*a;
%Derek - interpolate to gauss nodes
gRoota = gauss.interp*roota;


% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)'; 
for k1=1:K 
  if(~mod(k1,200)) k1, end;

  % Location of k1'th diagonal block entries in OP matrix
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

  % Extract local mass matrix and cubature weights
  locmm = cub.mm(:,:,k1);
  cw  = spdiags(cub.W (:,k1), 0, cub.Ncub, cub.Ncub);
  %Derek - define inverse of cubature weights matrix.
  cwinv = spdiags(1./cub.W(:,k1), 0, cub.Ncub, cub.Ncub);
  %Make diagonal matrix for local 'roota' 
  cRootaMat = spdiags(cRoota(:,k1), 0, cub.Ncub, cub.Ncub);
  cubAMat = spdiags(cubA(:,k1),0,cub.Ncub,cub.Ncub);
  
  % Evaluate derivatives of Lagrange basis functions at cubature nodes
  [cDx, cDy]  = PhysDmatrices2D(x(:,k1), y(:,k1), cub.V);

  % Evaluate local stiffness matrix
  %OP11 = cDx'*cw*cDx + cDy'*cw*cDy;
  % Derek:
  OP11 = cDx'*cw*cRootaMat*cwinv*cRootaMat*cw*cDx + ...
         cDy'*cw*cRootaMat*cwinv*cRootaMat*cw*cDy;
  %OP11 = cDx'*cw*cubAMat*cDx  + ...
  %       cDy'*cw*cubAMat*cDy;
  
  % Build element-to-element parts of stiffness matrix for element k1
  for f1=1:Nfaces

    % Find neighbor
    k2 = EToE(k1,f1); f2 = EToF(k1,f1);

    idsM = (f1-1)*NGauss+1:f1*NGauss;
    
    % Extract Lagrange basis function -> Gauss node interpolation matrix
    gVM = gauss.finterp(:,:,f1);
    gVP = gauss.finterp(:,:,f2);    gVP = gVP(NGauss:-1:1,:);

    % Evaluate spatial derivatives of  Lagrange basis function at Gauss nodes
    [gDxM, gDyM] = PhysDmatrices2D(x(:,k1), y(:,k1),gVM);
    [gDxP, gDyP] = PhysDmatrices2D(x(:,k2), y(:,k2),gVP);

    % Evaluate normals at Gauss nodes on face
    gnx = spdiags(gauss.nx(idsM, k1), 0, NGauss, NGauss);
    gny = spdiags(gauss.ny(idsM, k1), 0, NGauss, NGauss);
    gw  = spdiags(gauss.W(idsM, k1),  0, NGauss, NGauss);  %Gauss weights
    % Derek - make face matrix diag matrix for roota
    gRootaMat = spdiags(gRoota(idsM,k1), 0, NGauss, NGauss);
    %Derek - define inverse of Gaussian weights matrix.
    gwinv = spdiags(1./gauss.W(idsM,k1), 0, NGauss, NGauss);
    
    % Compute normal derivatives of Lagrange basis functions at Gauss nodes
    gDnM = gnx*gDxM + gny*gDyM;
    gDnP = gnx*gDxP + gny*gDyP;

    % Locate global numbers of Lagrange nodes in neighboring element 
    cols2 = ones(Np,1)*((k2-1)*Np+1:k2*Np); 

    % Find minimum height of two elements sharing this face
    hinv = max(Fscale( 1 + (f1-1)*Nfp, k1), Fscale( 1 + (f2-1)*Nfp, k2));    

    % Set penalty scaling
    %gtau = 20*(N+1)*(N+1)*hinv; %orig
    gtau = 2*100*100*(N+1)*(N+1)*hinv; %derek.
    %note: changing this '20' to '100' seems to make
    %the inhomog. dirichlet BC case fail. odd.

    
    %keyboard;
    % Determine type of face
    switch(BCType(k1,f1))
      case {Dirichlet}
        % Dirichlet boundary face variational terms
	%OP11 = OP11 + ( gVM'*gw*gtau*gVM - gVM'*gw*gDnM - gDnM'*gw*gVM); %orig
    
    OP11 = OP11 + (gtau*gVM'*gw*gRootaMat*gVM - gVM'*gw*gRootaMat*gDnM -gDnM'*gw*gRootaMat*gwinv*gRootaMat*gw*gVM);
    %Derek: I don't think we want to change the BC terms for sdd, 
    %may want to check this in nodal codes too.

      case {Neuman}
        % Do nothing
      case {Wall}
        % Do nothing, same as Neuman     
      otherwise
	% Interior face variational terms for stiffness matrix
	%OP11 = OP11 + 0.5*( gtau*gVM'*gw*gVM - gVM'*gw*gDnM - gDnM'*gw*gVM ); %orig.
    %Derek:
    OP11 = OP11 + 0.5*( gtau*gVM'*gw*gRootaMat*gVM - gVM'*gw*gRootaMat*gDnM - gDnM'*gw*gRootaMat*gwinv*gRootaMat*gw*gVM );
    
	%OP12 =      - 0.5*( gtau*gVM'*gw*gVP + gVM'*gw*gDnP - gDnM'*gw*gVP ); %orig.
    %Derek:
    OP12 =      - 0.5*( gtau*gVM'*gw*gRootaMat*gVP + gVM'*gw*gRootaMat*gDnP - gDnM'*gw*gRootaMat*gwinv*gRootaMat*gw*gVP );

        % Store self-neighbor interaction term in global stiffness matrix
        OP(entries(:), :) = [rows1(:), cols2(:), OP12(:)];
	entries = entries + Np*Np;
    end 
  end   

  % Store k1'th self-self interaction term in global stiffness matrix
  OP(entries(:), :)   = [rows1(:), cols1(:), OP11(:)];
  MM(entriesMM(:), :) = [rows1(:), cols1(:), locmm(:)];
  entries = entries + Np*Np; entriesMM = entriesMM + Np*Np;
end  

% Convert OP and MM from coordinate storage format to Matlab's intrinsic sparse matrix format
OP = OP(1:max(entries)-Np*Np,:); 
OP = myspconvert(OP, Np*K, Np*K, 1e-15); MM = myspconvert(MM, Np*K, Np*K, 1e-15);
return
