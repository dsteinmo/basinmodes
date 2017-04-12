%try to refine just the corners.
%this does the refinement in a stupid way
%and furthermore, doesn't seem to treat BC's properly.
%a work in progress.

function [Nv,VX, VY, K, EToV,BCType] = GenBox_refine()

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

fd = inline('drectangle(p,-1,1,-1,1)','p');
fh = @huniform;
h0=0.16;
%h0=0.1;
%h0=0.25;
%h0 = 0.16; %h0=.1 gives 800 elements, 0.25 -> 128, 0.2 ->200
           %h0=.16 for the 288 element box mesh
Bbox = [-1 -1; 1 1];
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

%refinement stuff.
    %axis([-1 -0.65 -1 -0.65]);
    axis([0.65 1 0.65 1]);
    axis on;
    refineflag =1;
    while refineflag == 1
        title('Pan to region that needs refinement. Press the ''any'' key when finished.');
        pause;
        figure(1);
        
        title('select a boundary vertex on a boundary element for refinement...');
        [myx,myy] = ginput(1);
        mytol=1e-2; %10m tolerance
        dists = sqrt((VX-myx).^2 + (VY-myy).^2);
        myV = find(dists<mytol);
        if length(myV) ~= 1
            title('No vertex selected, assuming you''re done refining...');
            refineflag=0;
        else
            title(['You clicked on vertex number: ' num2str(myV)]);
            myelement = mod(find(EToV==myV),K);
            myelement(myelement==0) = K;
            %if myelement == 0
            %    myelement = K;
            %end
            if length(myelement) ~= 1
                %title('Error: You picked a vertex with multiple elements. Click corner vertices only. (for now)');
                title('Vertex corresponds to multiple elments. Click the other boundary vertex for disambiguity');
                [myx2,myy2] = ginput(1);
                mytol=1e-2; 
                dists = sqrt((VX-myx2).^2 + (VY-myy2).^2);
                myV2 = find(dists<mytol);
                localEToV = EToV(myelement,:);
                numpossible = length(localEToV);
                theind = mod(find(myV2 == localEToV),numpossible);
                if theind == 0
                    theind = numpossible;
                end
                
                if length(theind) == 1
                      disp(['you selected element number: ' num2str(myelement(theind))]);
                      localEToV = EToV(myelement(theind),:);
                      hold on;
                     plot(VX([localEToV localEToV(1)]),VY([localEToV localEToV(1)]),'-y','linewidth',2);
                     hold off; 
                     
                      %... work in progress. this is an easy refine as
                      %no neighbours are effected
                      
                else
                      title('Error: Click on boundary vertices only!');
                      pause(2);
                end
                
                
                
            else
                disp(['Element number: ' num2str(myelement)]);    
                localEToV = EToV(myelement,:); %1x3 array of the element we want to refine
                myVind = find(localEToV==myV); %index of myV in localEToV
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

                title('done refining that element and its affected neighbour.');
                pause(2);

            end

        end 
    end %end refinement stuff while loop

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

[EToE, EToF] = tiConnect2D(EToV);% [myx,myy] = ginput(1);

%This splits a triangle into four by drawing a joining the midpoints
%of the old triangle. My guess is that it fails if two elements
%lie along an edge. Where are you now Engsig-Karup!?!

%figure(17);
%PlotMesh2D;

%allocate BCType table.
BCType = 0*EToE;

%Insert the correct BC codes for boundaries
%BCType = CorrectBCTable(EToV,BCType,nodesInner,Wall,K);

%Note: original "CorrectBCTable" will fail if more than
%1 triangle edge is on the boundary. so in this script, use _v2
%BCType = CorrectBCTable(EToV,BCType,nodesOuter,Neuman,K);
BCType = CorrectBCTable_v2(EToV,VX,VY,BCType,fd,Neuman);


return

