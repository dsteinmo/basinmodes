clear;
close all;
setPaths;
Globals2D;

%------------------------------------------------------------------- 
N = 2; % Polynomial order used for approximation
MESH_FILE = 'ontario_gmsh.mat';
NUM_H_REFINES = 0;

USEMEANDEPTH = false; %use mean depth instead of full bathymetry?
 
f = 2*(2*pi/(3600*24))*sin(43.7*pi/180); %43.7 deg. lattitude.
g = 9.81;

numpot = 100;  %number of potential basis functions
numstrm = 100; %number of streamfunction basis functions

DUMP_TO_FILE = false;
%-------------------------------------------------------------------

basinmodes_curved;