function BCType = CorrectBCTable_derek(EToV,VX,VY,node,edge,BCType,BCcode,hack)

%modified by steinmoeller to not require the global vars.
%and to work without a distance function

%this script can properly address elements with nodes on more than one wall,
%or two faces along a 'corner' boundary.

%Globals2D;

if ~exist('hack','var')
    hack = false;
end

%keyboard;
%VNUM = [1 2;2 3;3 1]; % face orientations
BCType = zeros(size(BCType));



pxc = 0.5*(VX(EToV)+VX(EToV(:,[2 3 1])));
pyc = 0.5*(VY(EToV)+VY(EToV(:,[2 3 1])));

midps = [pxc(:) pyc(:)];
tol = 1e-2;  %was 1e-4

%try this hack to fix colinearity bug with 'findedge'...
%put in dummy point:
if hack == true
    midps = [midps;1e9 1e9];
    edgenum = findedge(midps,node,edge,tol); %normal code.
    edgenum = edgenum(1:end-1);
else% drop dummy point's edge/flag. end hack.
    edgenum = findedge(midps,node,edge,tol);
end
%if the face's midpoint has an edgenumber, then it lies on the boundary!
idx = edgenum ~=0;

%dc  = abs(fd([pxc(:) pyc(:)])); % distances to boundaries from face centers
%tol = 1e-4; % tolerance
%idx = find(dc<tol);
%idx = dc<tol;

%BCType(idx) = BCcode;
BCType(edgenum>0) = BCcode;
return
