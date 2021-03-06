function [Nv,VX, VY, K, EToV,BCType] = GenDonut()

% function [VX, VY, K, EToV] = MeshGenDistMesh2D()
% Purpose  : Generate 2D square mesh using DistMesh;   
% By Allan P. Engsig-Karup

% Parameters to set/define
%    fd     Distance function for mesh boundary
%    fh     Weighting function for distributing elements
%    h0     Characteristic length of elements
%    Bbox   Bounding box for mesh
%    param  Parameters to be used in function call with DistMesh

fd = inline('drectangle(p,-1,1,-1,1)','p');
%fh = @huniform;
h0 = 0.1;
Bbox = [-1 -1; 1 1];
param = [];

% Call distmesh
% [Vert,EToV]=distmesh2d(fd,fh,h0,Bbox,param);  %this seems messed up, rectangle

% fd=inline('sqrt(sum(p.^2,2))-1','p');
% [Vert,EToV]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);  %unit circle

fd = inline('-0.3 +abs(0.7-sqrt(sum(p.^2,2)))');
% fh = inline('0 +sqrt(sum(p.^2,2)-0.1) ');
% fh = inline('0.05 + 0.15*sqrt(sum((p+1).^2,2))'); figure out how to use
% this crap, fplot?
fh = @huniform;

figure(1);
[Vert,EToV] = distmeshnd(fd,fh,0.1,[-1,-1;1,1],[]);

fd_inner = inline('sqrt(sum(p.^2,2))-0.4','p');
fd_outer = inline('sqrt(sum(p.^2,2))-1','p');

tol = 1e-8;
nodesInner = find(abs(fd_inner(Vert))<tol);
nodesOuter = find(abs(fd_outer(Vert))<tol);

%BC flags
In=1;
Out=2;
Wall=3;

VX = Vert(:,1); VY = Vert(:,2);
Nv = length(VX); K  = size(EToV,1);


% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);

% Build connectivity matrix
[EToE, EToF] = tiConnect2D(EToV);


%allocate BCType table.
BCType = 0*EToE;

%Insert the correct BC codes for boundaries
BCType = CorrectBCTable(EToV,BCType,nodesInner,Wall,K);

BCType = CorrectBCTable(EToV,BCType,nodesOuter,Wall,K);


%Need to do this to make vertex arrays consistent with main scripts.
VX = VX';
VY = VY';
return
