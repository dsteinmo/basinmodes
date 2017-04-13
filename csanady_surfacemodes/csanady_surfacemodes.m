clear;
close all;
currDir = pwd;
cd ../
root = pwd;
addpath(root);
setPaths;
cd(currDir);

Globals2D;

%------------------------------------------------------------------- 
N = 2; % Polynomial order used for approximation
MESH_FILE = '../circlemesh_nohole.mat';
NUM_H_REFINES = 0;

USEMEANDEPTH = false; %use mean depth instead of full bathymetry?
 
f = 2*(2*pi/(3600*24))*sin(43.7*pi/180); %43.7 deg. lattitude, (for model great lake, From Csanady 1967)
g = 9.81;

numpot = 150;  %number of potential basis functions
numstrm = 150; %number of streamfunction basis functions

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

% Analytic n=1,m=0 mode: -6.96 (Csanday 1967)
% Analytic n=1,m=0 mode:  7.80

% Num_H_REFINES = 0 (for next 3).
% N=1, numpot=150, numstrm=150. get: sigma_1/f = 6.9784, sigma_2/f = 7.8113 (very close to analytics)
% N=2, numpot=150, numstrm=150. get: sigma_1/f = 6.9671, sigma_2/f = 7.8035 (even better!)
% N=3, numpot=150, numstrm=150. get: sigma_1/f = 6.9668, sigma_2/f = 7.8038

%NUM_H_REFINES = 1
% N=2, numpot=150, numstrm=150. get: sigma_1/f = 6.9668, sigma_2/f = 7.8038