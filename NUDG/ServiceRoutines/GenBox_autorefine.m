%try to refine just the corners.
%this does the refinement in a stupid way
%and furthermore, doesn't seem to treat BC's properly.
%a work in progress.

function [Nv,VX, VY, K, EToV,BCType] = GenBox_autorefine()

%Globals2D;

% function [VX, VY, K, EToV] = MeshGenDistMesh2D()
% Purpose  : Generate 2D square mesh using DistMesh;   
% By Allan P. Engsig-Karup

% Parameters to set/define
%    fd     Distance function for mesh boundary
%    fh     Weighting function for distributing elements
%    h0     Characteristic length of elements
%    Bbox   Bounding box for mesh
%    param  Parameters to be used in function call with DistMesh

%fd = inline('drectangle(p,-1,1,-1,1)','p');
fd = inline('drectangle(p,0,1,0,1.01)','p');
fh = @huniform;
h0=0.05;
%h0=0.16;
%h0=0.1;
%h0=0.25;
%h0 = 0.16; %h0=.1 gives 800 elements, 0.25 -> 128, 0.2 ->200
           %h0=.16 for the 288 element box mesh
%Bbox = [-1 -1; 1 1];
Bbox = [0 0; 1 1.01];
param = [];


figure(1);
[Vert,EToV] = distmeshnd(fd,fh,h0,Bbox,param);

%fd_outer = inline('drectangle(p,-1,1,-1,1)','p');
fd_outer = inline('drectangle(p,0,1,0,1.01)','p');

tol = 1e-2;
%nodesInner = find(abs(fd_inner(Vert))<tol);
nodesOuter = find(abs(fd_outer(Vert))<tol);

%BC flags
In=1;
Out=2;
Wall=3;
Dirichlet=6;
Neuman=7;

VX = Vert(:,1); VY = Vert(:,2);
Nv = length(VX); K  = size(EToV,1);



% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);

%get first attempt at BC table.
BCType= zeros(K,3);
BCType = CorrectBCTable(EToV,BCType,nodesOuter,Neuman,K);
%find elements that have "3" faces on the boundary (means CorrectBCtable
%failed), then split.
badels=[];
for jj=1:K
    if BCType(jj,1) > 0 && BCType(jj,2) > 0 && BCType(jj,3) > 0
        badels = [badels; jj];
    end
end
disp('hi');
disp(length(badels));
disp(badels);

%find out which side is the hypoteneuse, and split accordingly.
for jj=1:length(badels)
    myelement = badels(jj);
    V1 = EToV(myelement,1);
    V2 = EToV(myelement,2);
    V3 = EToV(myelement,3);
    
    dists = [sqrt((VX(V1)-VX(V2)).^2 + (VY(V1)-VY(V2)).^2);
            sqrt((VX(V1)-VX(V3)).^2 + (VY(V1)-VY(V3)).^2);
            sqrt((VX(V2)-VX(V3)).^2 + (VY(V2)-VY(V3)).^2)];
    [aa bb] = max(dists);  %figure out which edge is hypoteneuse
    
    %the non-hypoteneuse vertex is the one on the boundary
    if bb == 1
        myVind = 3;
    elseif bb== 2
        myVind = 2;
    elseif bb== 3;
        myVind = 1;
    else
        disp('something bad happened');
        return;
    end
    
    localEToV = EToV(myelement,:);
    myV = localEToV(myVind);
    otherinds = find(localEToV ~= myV); %other two vertex indices
    other(1) = localEToV(otherinds(1)); %get global vertex numbers
    other(2) = localEToV(otherinds(2)); %of other two

    %calculate midpoint of other two vertices
    midpointx = 0.5*(VX(other(1)) + VX(other(2)));
    midpointy = 0.5*(VY(other(1)) + VY(other(2)));

    %add midpoint as new vertex to global vertex table
    VX=[VX;midpointx]; VY=[VY;midpointy];
    %midpoint's is at the end of the table, so its index is
    %length(table)
    midind = length(VX);

    %draw new line segment on mesh to indicate that refinement is
    %being done
    hold on;
    plot([midpointx VX(myV)],[midpointy VY(myV)],'-b'); drawnow;
    hold off;

    %find the element our 'split-face' is shared with and refine it
    %appropriately
    myinds = mod(find(EToV == other(1)),K);
    myinds(myinds==0) = K; %need this hack b/c of matlab.
    %if myinds == 0 
    %    myinds =K;
    %end

    myinds = myinds(myinds ~= myelement); %remove the element we started with from the contenders
    numpossible = length(myinds);
    possibleelements = EToV(myinds,:);
    theind = mod(find(possibleelements == other(2)),numpossible);
    if theind == 0 %need this hack to get around matlab no zero index problems
        theind = numpossible;
    end
    otherelement = myinds(theind);
    othervertices = EToV (otherelement,:);       %global vertex numbers of the 'other element'
    %othervertices = possibleelements(theind,:); %global vertex numbers of the 'other element'

    tmpinds = find(othervertices ~= other(2));
    newmyVind = othervertices(tmpinds) ~= other(1);
    newmyVind = tmpinds(newmyVind);

    newmyV = othervertices(newmyVind);

    %draw new line segment on mesh to indicate that refinement is
    %being done
    hold on;
    plot([midpointx VX(newmyV)],[midpointy VY(newmyV)],'-b'); drawnow;
    hold off;

    %modify old element indices to include midpoint
    EToV(myelement,:) = [midind other(1) myV];
    %add a new element to the connectivity table
    newel = [midind myV other(2)];
    EToV = [EToV; newel];

    %modify old element indices to include midpoint
    EToV(otherelement,:) = [newmyV other(1) midind];
    %add a new element to the connectivity table
    newel = [newmyV midind other(2)];
    EToV = [EToV; newel];


    %increment number of vertices and increment twice the num. of elements
    Nv=Nv+1
    K=K+2

    %highlight our refined elemenets (repeat first vertex to close
    %off triangles)
    hold on;
    plot(VX([EToV(myelement,:) EToV(myelement,1)]),VY([EToV(myelement,:) EToV(myelement,1)]),'-r','linewidth',2);
    plot(VX([EToV(otherelement,:) EToV(otherelement,1)]),VY([EToV(otherelement,:) EToV(otherelement,1)]),'-g','linewidth',2);
    plot(VX([EToV(end-1,:) EToV(end-1,1)]),VY([EToV(end-1,:) EToV(end-1,1)]),'-m','linewidth',2);
    plot(VX([EToV(end,:) EToV(end,1)]),VY([EToV(end,:) EToV(end,1)]),'-c','linewidth',2);
    drawnow;
    hold off;

    
end


% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);
    

%Need to do this to make vertex arrays consistent with main scripts.
VX = VX';
VY = VY';

%[EToE, EToF] = tiConnect2D(EToV);% [myx,myy] = ginput(1);

%This splits a triangle into four by drawing a joining the midpoints
%of the old triangle. My guess is that it fails if two elements
%lie along an edge. Where are you now Engsig-Karup!?!

%figure(17);
%PlotMesh2D;

%allocate BCType table.
%BCType = 0*EToE;

%Insert the correct BC codes for boundaries
%BCType = CorrectBCTable(EToV,BCType,nodesInner,Wall,K);

%Note: original "CorrectBCTable" will fail if more than
%1 triangle edge is on the boundary. so in this script, use _v2
BCType = zeros(K,3);
BCType = CorrectBCTable(EToV,BCType,nodesOuter,Neuman,K);
%BCType = CorrectBCTable_v2(EToV,VX,VY,BCType,fd,Neuman);

hold on;
plot(VX(nodesOuter),VY(nodesOuter),'.r');
drawnow;
hold off;


return

