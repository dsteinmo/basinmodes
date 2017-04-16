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
ANALYTIC_DEPTH = true;

f = 2*(2*pi/(3600*24))*sin(43.7*pi/180); %43.7 deg. lattitude, (for model great lake, From Csanady 1967)
g = 9.81;

numpot = 500;  %number of potential basis functions
numstrm = 500; %number of streamfunction basis functions

DUMP_TO_FILE = false;

Hepi = 15;       %Csanady calls "h",  epilimnion  thickness (m)
Hhyp = 60;       %Csanady calls "h'", hypolimnion thickness (m)

deltarho = 1.74; % (converted to kg/m^3)
rhoprime =  1.00e+03; %hypolimnion density

epsilon = deltarho/rhoprime; %relative density difference
gprime = epsilon*g;  %reduced gravity

H1 = Hepi+Hhyp;  %total lake depth, (equivalent depth for barotropic mode)
H2 = epsilon*((Hepi*Hhyp)/(Hepi+Hhyp)); %equivalent depth for internal mode.

% Equivalent Depth is a constant. We give a trivial function-handle to make the main code happy.
H_analytic = @(x,y)H2 + 0*x;

CURVED_IMPLEMENTATION = true;
%-------------------------------------------------------------------

if CURVED_IMPLEMENTATION == true
    basinmodes_curved;
else
    basinmodes;
end
%Baroclinic Modes (sigma/f). DG N=2, numpot=500, numstrm=500, hrefines=0
%--------------
%Gravest Poincare mode

%n   m     Analytic (Csanady1967)         DG
%-   -     ----------------------       -------
%1   0            1.013                  1.012667

%Kelvin Modes (sigma/f)
%------------

% n        Analytic (Csanady1967)         DG (hrefines=0)			DG (hrefines=1)
%---       ----------------------       ---------					---------------
% 1               -0.0689               -0.069239					-0.069256
% 2               -0.1373               -0.138464					-0.138278
% 4               -0.2747               -0.276887					-0.276629
% 9               -0.6211               -0.62270					-0.621931
% 14              -0.9333               -0.96812					-0.966997 *(~3.6% error!)