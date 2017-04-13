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

c1 = sqrt(g*H1);
c2 = sqrt(g*H2);

r0 = 6.75e4;
%Lamb's Beta parameter
beta=40;
f = sqrt(beta)*(c2/r0);

H_analytic = @(x,y)H2*(1-(x.^2 + y.^2)/r0^2) + 0.005*H2;

basinmodes_curved;

% beta: 40
% --------
% numstrm+82  (582) => (s=1, n=3, cw),  1.0482
% numstrm-145 (355) => (s=1, n=3, ccw), 0.048223

% numstrm+101 (601) => (s=1, n=5, cw),  1.1865
% numstrm-183 (317) => (s=1, n=5, ccw), 0.037690
% numstrm+93  (593) => (s=1, n=5, ccw), 1.1503
