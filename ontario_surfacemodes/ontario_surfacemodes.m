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
MESH_FILE = 'ontario_gmsh2.mat';
NUM_H_REFINES = 0;

USEMEANDEPTH = false; %use mean depth instead of full bathymetry?
ANALYTIC_DEPTH = false;
 
f = 2*(2*pi/(3600*24))*sin(43.7*pi/180); %43.7 deg. lattitude.
g = 9.81;

numpot = 100;  %number of potential basis functions
numstrm = 100; %number of streamfunction basis functions

DUMP_TO_FILE = false;

CURVED_IMPLEMENTATION = true;
%-------------------------------------------------------------------

if CURVED_IMPLEMENTATION == true
    basinmodes_curved;
else
    basinmodes;
end