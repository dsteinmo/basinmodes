%asc2mesh


%read parameters from header file
fid = fopen('newont-6538.asc','r');
%had to get this by plotting the bathymetry, post-projection
%R_E = 6371000; %earth radius in m;

%c = fread(fid,12,'uint8=>char');
l = fgetl(fid); Nlong = str2double(l(8:end));
l = fgetl(fid); Nlat = str2double(l(8:end));
l = fgetl(fid); longstart = str2double(l(12:end));
l = fgetl(fid); latstart = str2double(l(12:end));
l = fgetl(fid); dtheta = str2double(l(10:end)); %grid spacing in degrees
l = fgetl(fid); nodata = str2double(l(15:end));
fclose(fid);
% %s = fread(fid,120,'uint8=>char',4)
% s = fread(fid,9,'uint8=>char')
% s= strcat(s')
% str2num(s)

%load bathymetry from a version of the file with header deleted.
%(Note: could try doing this with single file & fread);
if ~exist('bathy','var')
    disp('loading bathymetry data...');
    bathy = load('newont.asc');
    disp('done');
end


[Nlat Nlong] = size(bathy);

%build longitude/latitude grid from this data
long = longstart + (0:Nlong-1)*dtheta;
%lat = latstart + (0:Nlat-1)*dtheta; %had this
lat = latstart - (0:Nlat-1)*dtheta;

[long2,lat2] = meshgrid(long,lat);

%mstruct = defaultm('tranmerc');
%mstruct = defaultm('lambertstd');
 %mstruct = defaultm('eqdconicstd'); %makes lake too small.
%mstruct = defaultm('cassinistd');
mstruct = defaultm('polyconstd'); %schwab, basically same as tranverse mercator
%mstruct = defaultm('tranmerc'); %that noaa website used transverse mercator
%mstruct = defaultm('eqdazim');
%mstruct = defaultm('mercator'); %makes the lake too big.
mstruct.geoid = almanac('earth','grs80','meters');

longmin = min(long); longmax = max(long);
latmin = min(lat); latmax = max(lat);

mstruct.maplonlimit = [longmin longmax];
mstruct.maplatlimit = [latmin latmax];
mstruct = defaultm(mstruct);


%get (x,y)-coords for bounding box.
[xmin ymin] = projfwd(mstruct,latmin,longmin);
[xmax ymax] = projfwd(mstruct,latmax,longmax);

%make a 1024x1024 cartesian grid based on these bounds
Ny= Nlat; Nx = Nlong;

Nxout = 1024; dx = (xmax-xmin)/Nxout;
Nyout = 1024; dy = (ymax-ymin)/Nyout;

xc = xmin+ (0:Nxout-1)*dx;
yc = ymin+ (0:Nyout-1)*dy;

[xx,yy] = meshgrid(xc,yc);

%replace any positive land values, 'nodata values', or nans with 0 for
%land.
bathy(isnan(bathy) | bathy == nodata | bathy  > 0) = 0;

%find what long/lat values our cartesian grid corresponds to.
[mylat mylong] = projinv(mstruct,xx,yy);
%interpolate from the values of the input data to this set of long/lats
%hence giving us the data at the points of our cartesian grid.
extrapval = 0;
bathycart = interp2(long2,lat2,bathy,mylong,mylat,'linear',extrapval);

%now shift coordinates so bottom left is (0,0)
xx = xx - min(xx(:));
yy = yy - min(yy(:));

% Now, smooth the bathymetry
bzsmooth = bathycart;

for k=1:8 %8
    bzsmooth = conv2(bzsmooth,[1/8 3/4 1/8]'*[1/8 3/4 1/8],'same');
end

land = double(bzsmooth > -0.5); % Find the 'land' cells
for k=1:3 %3
    % Smooth that into the water cells, to find near-land cells.
    land = conv2(land,ones(3),'same');
end
% Invert that map, to find water that is -not- disconnected
water = double(land == 0);
for k=1:3  %3
    water = conv2(water,ones(3),'same');
end

bzsmooth((water == 0)) = 0;
bzsmooth(bzsmooth > -0.5) = 0;


%to help with meshing restrict to coarser mesh

Nxmesh = 128; dxm = max(xx(:))/Nxmesh;
Nymesh = 128; dym = max(yy(:))/Nymesh;
xm = (0:Nxmesh-1)*dxm;
ym = (0:Nymesh-1)*dym;

[xxm,yym] = meshgrid(xm,ym);

bzmesh = interp2(xx,yy,bathycart,xxm,yym);

%bzmesh = bathycart;

for k=1:8 %8
    bzmesh = conv2(bzmesh,[1/8 3/4 1/8]'*[1/8 3/4 1/8],'same');
end

land = double(bzmesh > -0.5); % Find the 'land' cells
for k=1:3 %3
    % Smooth that into the water cells, to find near-land cells.
    land = conv2(land,ones(3),'same');
end
% % Invert that map, to find water that is -not- disconnected
water = double(land == 0);
for k=1:3  %3
    water = conv2(water,ones(3),'same');
end

bzmesh((water == 0)) = 0;
bzmesh(bzmesh > -0.5) = 0;

%done smoothing bathymetry, now try generating node/edge connectivity

% generate the zero contour
hold off;
figure(16);
[C,h] = contour(xxm,yym,bzmesh,[0 0],'-','linewidth',2);
%keyboard;

%return;


% check that we actually made contour lines
if size(C,2) > 1
    stillContours = true;
else
    stillContours = false;
    disp('warning: no contour lines generated');
end

node = [];
edge = [];
currStart = 1;
oldNumNodes = 0;

% generate node and node connectivity data
while stillContours
    
    % number of nodes in current contour line
    currNumNodes = C(2,currStart);
    
    % x- and y-coordinates of nodes
    currX = C(1,currStart+(1:currNumNodes))';
    currY = C(2,currStart+(1:currNumNodes))';
    
    % add nodes to node list
    node = [node; currX, currY];
    
    % get the new number of nodes
    newNumNodes = size(node,1);
    
    % generate edge connectivity data
    edge = [edge; 
            ((oldNumNodes+1):(newNumNodes-1))' ((oldNumNodes+2):newNumNodes)';
            newNumNodes oldNumNodes+1];
    
    % compute new starting position
    currStart = currStart + currNumNodes + 1;
    
    oldNumNodes = newNumNodes;
    
    % check if we've reached the end yet
    if currStart >= size(C,2)
        stillContours = false;
    end
end %end while


disp('drawing node/edge connectivity...');

%keyboard;
figure(9); clf;

hold on;
for jj=1:length(edge)
    plot([node(edge(jj,1),1) node(edge(jj,2),1)],[node(edge(jj,1),2) node(edge(jj,2),2)],'-k');
end
hold off;
disp('done');

[Vert,EToV] = mesh2d(node,edge,[]);
axis on;
grid on; 
xlabel('x (m)'); ylabel('y (m)');
axis tight;
dimp = size(Vert); dimt = size(EToV);
Nv = dimp(1); K = dimt(1);

%stuff below for DG
VX = Vert(:,1); VY = Vert(:,2);

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);
%done reordering

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

%calculate # of DoF's if we were asked to
if exist('N','var')
   dof = nchoosek(N+2,2)*K;
   disp(['DG simulation with order N=' num2str(N) ' basis polynomials will have ' num2str(dof) ' degrees of freedom.']);
end



%find all boundary nodes (nodesOuter). Ain't nothin to it, but to do it.
tol = 1e-8;
edgenum = findedge(Vert,node,edge,tol);
kk=1;
nodesOuter = [];
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

%allocate BCType table.
BCType = 0*EToE;

%Insert the correct BC codes for boundaries
Wall=3;
BCType = CorrectBCTable_derek(EToV,VX,VY,node,edge,BCType,Wall);

%Need to do this to make vertex arrays consistent with main scripts.
VX = VX';
VY = VY';

%shift things so that southwest is (0,0)
xmin = min(VX(:));
ymin = min(VY(:));

VX = VX - xmin;
VY = VY - ymin;
xx = xx - xmin;
yy = yy - ymin;

figure(19);
pcolor(xx,yy,bzsmooth); shading flat; colorbar;

save('lake_ontario_mesh_coarse.mat','Nv','VX','VY','K','EToV','BCType','xx','yy','bzsmooth');
