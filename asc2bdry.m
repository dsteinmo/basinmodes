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


% figure(8);
% pcolor(xx,yy,bzsmooth); shading flat; colorbar; 


%to help with meshing restrict to coarser mesh

Nxmesh = 150; dxm = max(xx(:))/Nxmesh;
Nymesh = 150; dym = max(yy(:))/Nymesh;
xm = (0:Nxmesh-1)*dxm;
ym = (0:Nymesh-1)*dym;

[xxm,yym] = meshgrid(xm,ym);

bzmesh = interp2(xx,yy,bathycart,xxm,yym);

%bzmesh = bathycart;

for k=1:1 %8
    bzmesh = conv2(bzmesh,[1/8 3/4 1/8]'*[1/8 3/4 1/8],'same');
end

land = double(bzmesh > -0.5); % Find the 'land' cells
for k=1:2 %3
    % Smooth that into the water cells, to find near-land cells.
    land = conv2(land,ones(3),'same');
end
% % Invert that map, to find water that is -not- disconnected
water = double(land == 0);
for k=1:2  %3
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
conts=[];
n=1; %number of contours

edgearray = cell(1,1);

% generate node and node connectivity data
while stillContours
    
    % number of nodes in current contour line
    currNumNodes = C(2,currStart);
    
    % x- and y-coordinates of nodes
    currX = C(1,currStart+(1:currNumNodes))';
    currY = C(2,currStart+(1:currNumNodes))';
    
    %Derek: store individual contours in struct array for later splining
    conts(n).x = currX;
    conts(n).y = currY;
    
    
    % add nodes to node list
    node = [node; currX, currY];
    
    % get the new number of nodes
    newNumNodes = size(node,1);
    
    % generate edge connectivity data
    edge = [edge; 
            ((oldNumNodes+1):(newNumNodes-1))' ((oldNumNodes+2):newNumNodes)';
            newNumNodes oldNumNodes+1];
        
    edgearray{n} = [((oldNumNodes+1):(newNumNodes-1))' ((oldNumNodes+2):newNumNodes)';
            newNumNodes oldNumNodes+1];
    
    % compute new starting position
    currStart = currStart + currNumNodes + 1;
    
    oldNumNodes = newNumNodes;
    
    % check if we've reached the end yet
    if currStart >= size(C,2)
        stillContours = false;
    end
    n=n+1;
end %end while


%generate spline interpolants based on contour data
bdryspline=[];
for j=1:length(conts)
    [ bdryspline(j).xs, bdryspline(j).ys, bdryspline(j).ts ] = ParametricSpline(conts(j).x , conts(j).y);
end


depthdata.H = -bzsmooth; %convert to positive values for shallow-water solvers.
mindepth = 1;
%depthdata.H(depthdata.H<abs(mindepth))=abs(mindepth);
%for k=1:2
%    depthdata.H = conv2(depthdata.H,[1/8 3/4 1/8]'*[1/8 3/4 1/8],'same');
%end
depthdata.H(depthdata.H<abs(mindepth))=abs(mindepth);
depthdata.x = xx;
depthdata.y = yy;


disp('drawing node/edge connectivity...');

%keyboard;
figure(9); clf;

hold on;
for jj=1:length(edge)
    
    plot([node(edge(jj,1),1) node(edge(jj,2),1)],[node(edge(jj,1),2) node(edge(jj,2),2)],'-k');
    %jj
end
hold off;
disp('done');

