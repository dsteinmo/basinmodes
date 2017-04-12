function [VX, VY, K, EToV] = MeshGenDistMesh2D()

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
fh = @huniform;
h0 = 0.25;
Bbox = [-1 -1; 1 1];
param = [];

% Call distmesh
% [Vert,EToV]=distmesh2d(fd,fh,h0,Bbox,param);  %this seems messed up, rectangle

% fd=inline('sqrt(sum(p.^2,2))-1','p');
% [Vert,EToV]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);  %unit circle

fd = inline('-0.3 +abs(0.7-sqrt(sum(p.^2,2)))');
[Vert,EToV] = distmeshnd(fd,@huniform,0.1,[-1,-1;1,1],[]);

fd_inner = inline('sqrt(sum(p.^2,2))-0.4','p');
fd_outer = inline('sqrt(sum(p.^2,2))-1','p');

VX = Vert(:,1); VY = Vert(:,2);
Nv = length(VX); K  = size(EToV,1);

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);
return
