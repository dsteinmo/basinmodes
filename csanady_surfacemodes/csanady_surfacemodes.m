clear;
close all;
addpath '../';
setPaths;
Globals2D;

%------------------------------------------------------------------- 
N = 2; % Polynomial order used for approximation
MESH_FILE = 'circlemesh_nohole.mat'; 
NUM_H_REFINES = 0;

USEMEANDEPTH = false; %use mean depth instead of full bathymetry?
 
f = 2*(2*pi/(3600*24))*sin(43.7*pi/180); %43.7 deg. lattitude, (for model great lake, From Csanady 1967)
g = 9.81;

numpot = 100;  %number of potential basis functions
numstrm = 100; %number of streamfunction basis functions

DUMP_TO_FILE = false;

% We don't have any cartesian bathymetry data for idealized lake, so we fake it here.
depthdata = [];
Nx = 100;
Ny = 100;
r0 = 67.5e3;
xc = linspace(-r0, r0, Nx);
yc = linspace(-r0, r0, Ny);
[depthdata.x, depthdata.y] = meshgrid(xc, yc);
depthdata.H = 75.0*ones(Ny, Nx);
%-------------------------------------------------------------------

basinmodes_curved;
