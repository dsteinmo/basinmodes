
%filename should point to a mat-file that contains
%Nv, VX, VY, K, EToV, BCType
%as well as boundary node data 'node' and node-connectivity array "edge" to
%correct the BC-table after refining is done.

function [Nv,VX, VY, K, EToV,BCType] = myrefinetool(filename)
Globals2D; %stupid

mytol=1e-2; %click tolerance. adjust depending on length scales involved.

%load in mesh data. should have Nx,Vx,Vy,K and EToV in it.
load(filename);


%Plot the mesh with Warburton's crap

N=1;
StartUp2D; BuildMaps2D;
hold off; PlotMesh2D; axis tight;

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);


%refinement stuff.
    %axis([-1 -0.65 -1 -0.65]);
    %daxis([0.65 1 0.65 1]);
    axis on;
    refineflag =1;
    while refineflag == 1
        title('Zoom/Pan to region that needs refinement. Press the ''any'' key when finished.');
        pause;
        figure(1);
        
        title('select a boundary vertex on a boundary element for refinement...');
        [myx,myy] = ginput(1);
        dists = sqrt((VX-myx).^2 + (VY-myy).^2);
        myV = find(dists<mytol);
        if length(myV) ~= 1
            title('No vertex selected. Are you done refining? (Mouse button = no, Key=yes)');
            disp('No vertex selected. Are you done refining? (Mouse button = no, Key=yes)...');
            keydown = waitforbuttonpress;
            if keydown ~=0
                refineflag=0;
            end
        else
            hold on; plot(VX(myV),VY(myV),'.r');
            title(['You clicked on vertex number: ' num2str(myV)]);
            myelement = mod(find(EToV==myV),K);
            myelement(myelement==0) = K;
            if length(myelement) ~= 1 
                %title('Error: You picked a vertex with multiple elements. Click corner vertices only. (for now)');
                title('Vertex corresponds to multiple elments. Click the other boundary vertex for disambiguity');
                [myx2,myy2] = ginput(1);
                
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
                     
                     
                   
                     
%                      hold on;
%                      plot(VX([localEToV localEToV(1)]),VY([localEToV localEToV(1)]),'-y','linewidth',2);
%                      hold off; 
                     
                     %reshuffle things so we can apply old code.
                     other(1) = myV;
                     other(2) = myV2;
                     
                    myV = localEToV(localEToV ~= other(1) & localEToV ~= other(2));
                     
                     %calculate midpoint of other two vertices
                    midpointx = 0.5*(VX(other(1)) + VX(other(2)));
                    midpointy = 0.5*(VY(other(1)) + VY(other(2)));

                    
                    
                    %add midpoint as new vertex to global vertex table
                    VX=[VX midpointx]; VY=[VY midpointy];
                    %midpoint's is at the end of the table, so its index is
                    %length(table)
                    midind = length(VX);

                    %draw new line segment on mesh to indicate that refinement is
                    %being done
                    hold on;
                    plot([midpointx VX(myV)],[midpointy VY(myV)],'-k'); drawnow;
                    hold off;
                
                    %modify old element indices to include midpoint
                    EToV(myelement(theind),:) = [midind other(1) myV];
                    %add a new element to the connectivity table
                    newel = [midind myV other(2)];
                    EToV = [EToV; newel];
                    disp('refinement done...');
                    
                    K=K+1
                    Nv=Nv+1

                                    %highlight our refined elemenets (repeat first vertex to close
                %off triangles)
                hold on;
                plot(VX([EToV(myelement(theind),:) EToV(myelement(theind),1)]),VY([EToV(myelement(theind),:) EToV(myelement(theind),1)]),'-r','linewidth',2);
                plot(VX([EToV(end,:) EToV(end,1)]),VY([EToV(end,:) EToV(end,1)]),'-c','linewidth',2);
                drawnow;
                hold off;
                
                %this bit seems good.
                      
                elseif length(theind)==2
                      hold on; plot(VX(myV2),VY(myV2),'.r'); hold off;
                      title('Click the final vertex on the element for refinement...');
                      [myx3,myy3] = ginput(1);
                      dists = sqrt((VX-myx3).^2 + (VY-myy3).^2);
                      myV3 = find(dists<mytol);
                      
                      
                      localEToV = localEToV(theind,:);
                      foo = size(localEToV);
                      numpossible = foo(1);
                      theind2 = mod(find(myV3 == localEToV),numpossible);
                      if theind2 == 0
                        theind2 = numpossible;
                      end
                      
                      if length(theind2) ~=1
                          title('something bad happened.');
                          pause(2);
                      else
                          %disp(['You selected element number: ' num2str(myelement(theind(theind2)))]);
                          
                          localEToV = localEToV(theind2,:);
                          
                          %highlight the element
                          hold on;
                          plot(VX([localEToV localEToV(1)]),VY([localEToV localEToV(1)]),'-y','linewidth',2);
                          hold off; 
                          
                          title(['Element #: ' num2str(myelement(theind(theind2))) '. Click the vertex opposite the longest edge.']);
                          
                          [myx4,myy4] = ginput(1);
                          dists = sqrt((VX-myx4).^2 + (VY-myy4).^2);
                          myV4 = find(dists<mytol);
                          
                          
                          if length(myV4) ~= 1
                              title('Please click a vertex on the highlighted element.');
                              hold on;
                                plot(VX([localEToV localEToV(1)]),VY([localEToV localEToV(1)]),'-k','linewidth',2);
                              hold off; 
                              pause(2);
                          else
                          
                              %myV = find(myV4 == localEToV);
                              myV = localEToV(myV4 == localEToV);
                          
                              if length(myV) ~= 1
                                   title('Please click a vertex on the highlighted element.');
                                   hold on;
                                   plot(VX([localEToV localEToV(1)]),VY([localEToV localEToV(1)]),'-k','linewidth',2);
                                   hold off; 
                                   pause(2);
                              else  %good case.
                                    otherinds = find(localEToV ~= myV); %other two vertex indices
                                    other(1) = localEToV(otherinds(1)); %get global vertex numbers
                                    other(2) = localEToV(otherinds(2)); %of other two
                                    
                                    %keyboard;

                                    %calculate midpoint of other two vertices
                                    midpointx = 0.5*(VX(other(1)) + VX(other(2)));
                                    midpointy = 0.5*(VY(other(1)) + VY(other(2)));

                                    %add midpoint as new vertex to global vertex table
                                    VX=[VX midpointx]; VY=[VY midpointy];
                                    %midpoint's is at the end of the table, so its index is
                                    %length(table)
                                    midind = length(VX);

                                    %draw new line segment on mesh to indicate that refinement is
                                    %being done
                                    hold on;
                                    plot([midpointx VX(myV)],[midpointy VY(myV)],'-k'); drawnow;
                                    hold off;

                                    %find the element our 'split-face' is shared with and refine it
                                    %appropriately
                                    myinds = mod(find(EToV == other(1)),K);
                                    myinds(myinds==0) = K; %need this hack b/c of matlab.
                                    %if myinds == 0 
                                    %    myinds =K;
                                    %end

                                    myelement = myelement(theind(theind2));
                                    
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
                                    plot([midpointx VX(newmyV)],[midpointy VY(newmyV)],'-k'); drawnow;
                                    hold off;

                                    %keyboard;
                                    
                                    
                                    
                                    
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
                          
                          
                          
                      end
                      
                      
                else
                     title('Click on vertices that are shared by elements please!');                     
                end
                
                
                
            else  %we must be in the corner-vertex case
                disp(['Element number: ' num2str(myelement)]);    
                localEToV = EToV(myelement,:); %1x3 array of the element we want to refine
                %myVind = find(localEToV==myV); %index of myV in localEToV
                otherinds = find(localEToV ~= myV); %other two vertex indices
                other(1) = localEToV(otherinds(1)); %get global vertex numbers
                other(2) = localEToV(otherinds(2)); %of other two

                %calculate midpoint of other two vertices
                midpointx = 0.5*(VX(other(1)) + VX(other(2)));
                midpointy = 0.5*(VY(other(1)) + VY(other(2)));

                %add midpoint as new vertex to global vertex table
                VX=[VX midpointx]; VY=[VY midpointy];
                %midpoint's is at the end of the table, so its index is
                %length(table)
                midind = length(VX);

                %draw new line segment on mesh to indicate that refinement is
                %being done
                hold on;
                plot([midpointx VX(myV)],[midpointy VY(myV)],'-k'); drawnow;
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
                plot([midpointx VX(newmyV)],[midpointy VY(newmyV)],'-k'); drawnow;
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
%VX = VX';
%VY = VY';

[EToE, EToF] = tiConnect2D(EToV);% [myx,myy] = ginput(1);


%This splits a triangle into four by drawing a joining the midpoints
%of the old triangle. My guess is that it fails if two elements
%lie along an edge. Where are you now Engsig-Karup!?!

Vert = [VX' VY'];

%find all boundary nodes (nodesOuter). Ain't nothin to it, but to do it.
tol = 1e-8;
edgenum = findedge(Vert,node,edge,tol);
kk=1;
for jj = 1:length(edgenum)
    %if edgenum of point jj is nonzero, then it lies on the boundary,
    %so put it in list of boundary pts.
    if edgenum(jj) ~= 0  
        nodesOuter(kk) = jj;
        kk=kk+1;
    end        
end

hold on;
plot(VX(nodesOuter),VY(nodesOuter),'.g');
drawnow;
hold off;



%allocate BCType table.
BCType = 0*EToE;

disp('BC table not built yet. Please code this Derek.');

%Insert the correct BC codes for boundaries
%BCType = CorrectBCTable(EToV,BCType,nodesInner,Wall,K);

%Note: original "CorrectBCTable" will fail if more than
%1 triangle edge is on the boundary. so use _v2 if that happens.
BCType = CorrectBCTable(EToV,BCType,nodesOuter,Neuman,K);



return

