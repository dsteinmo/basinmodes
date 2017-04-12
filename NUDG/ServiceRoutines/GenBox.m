function [Nv,VX, VY, K, EToV,BCType] = GenBox()

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
%h0=0.1;
%h0 = 0.25;
h0=.16;
%h0 = 0.16; %h0=.1 gives 800 elements, 0.25 -> 128, 0.2 ->200
           %h0=.16 for the 288 element box mesh
Bbox = [0 0; 1 1];
param = [];




figure(1);
[Vert,EToV] = distmeshnd(fd,fh,h0,Bbox,param);

fd_outer = inline('drectangle(p,-1,1,-1,1)','p');

tol = 1e-2;
%nodesInner = find(abs(fd_inner(Vert))<tol);
%nodesOuter = find(abs(fd_outer(Vert))<tol);



%BC flags
In=1;
Out=2;
Wall=3;
Dirichlet=6;
Neuman=7;

VX = Vert(:,1); VY = Vert(:,2);
Nv = length(VX); K  = size(EToV,1);

%hold on;
%plot(VX(nodesOuter),VY(nodesOuter), '.r');
%hold off;


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
%BCType = CorrectBCTable(EToV,BCType,nodesInner,Wall,K);

%Note: original "CorrectBCTable" will fail if more than
%1 triangle edge is on the boundary. so in this script, use _v2
%BCType = CorrectBCTable(EToV,BCType,nodesOuter,Neuman,K);
BCType = CorrectBCTable_v2(EToV,VX,VY,BCType,fd,Neuman);

%Need to do this to make vertex arrays consistent with main scripts.
VX = VX';
VY = VY';
return

