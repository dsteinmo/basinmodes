clear;
close all;

node = [12 0;
        12 20;
        11 22;
        10 24;  %bay
        0 24;   %bay
        0 26;   %bay
        10 26;  
        11 28;
        12 30;
        30 30;
        30 0]*1e3;

edge = [1 2;
        2 3;
        3 4;
        4 5;
        5 6;
        6 7;
        7 8;
        8 9;
        9 10;
        10 11;
        11 1];

%make resolution finer in the bay
hdata.edgeh = [4 3e2; 
               5 3e2;
               6 3e2]; % Element size (300 m) on specified edges
                            
%highest resolution elsewhere (1.2 km)
hdata.hmax = 1.2e3;
[Vert,EToV] = mesh2d(node,edge,hdata);
axis on;
grid on; 
xlabel('x (m)'); ylabel('y (m)');
axis tight;
dimp = size(Vert); dimt = size(EToV);
Nv = dimp(1); K = dimt(1);
disp(['Mesh generated. ' num2str(Nv) ' vertices on ' num2str(K) ' elements.']);

%stuff below for DG
VX = Vert(:,1); VY = Vert(:,2);

%keyboard;

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);
%done reordering

% Build connectivity matrix
[EToE, EToF] = tiConnect2D(EToV);


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
plot(VX(nodesOuter),VY(nodesOuter),'.r');
drawnow;
hold off;

%Insert the correct BC codes for boundaries
Wall=3;
BCType = zeros(K,3);
BCType = CorrectBCTable(EToV,BCType,nodesOuter,Wall,K);

%Need to do this to make vertex arrays consistent with main scripts.
VX = VX';
VY = VY';

save('embayment2.mat','Nv','VX','VY','K','EToV','BCType');
    
    
