% Driver script for solving the 1D Poisson equation
Globals1D;

% Polynomial order used for approximation 
N =8;

% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(0,2*pi,10);


% Initialize solver and construct grid and metric
StartUp1D;

% Set RHS (Mass matrix times forcing)
%force = -(sin(2*x) + sin(x));
%f = J.*((invV'*invV)*force);  %to compare with exact solution

%what if RHS = d/dx of stuff ?
%stuff = cos(2*x)/2 + cos(x);
%dstuff = zeros(Nfp*Nfaces,K); dstuff(:) = (stuff(vmapM)-stuff(vmapP))/2.0;

%f = J.*((invV'*invV)*(rx.*(Dr*stuff) - LIFT*(Fscale.*dstuff) )); %Agrees with exact solution!

stuff = -cos(x).*(2*sin(x)+3);
f = J.*((invV'*invV)*stuff);

%set spatially dependent diffusivity
a = 2+sin(x);
ax = rx.*(Dr*a);

% Set up operator
%[A] = PoissonCstab1D();
[A] = HelmCstab1Dsdd(a);

% Solve Problem
solvec = A\f(:);
u = reshape(solvec,Np,K);

uexact = cos(x);

figure(1);
plot(x(:),u(:),x(:),uexact(:));

norm(u(:)-uexact(:),2)

%Now, try inhomogenous neuman bc with lifting operator
xR=pi/2;
xL=-pi/2;

a=1; b=-1; %neuman values

[Nv, VX, K, EToV] = MeshGen1D(xL,xR,100);

N=3;
% Initialize solver and construct grid and metric
StartUp1D;

%set spatially dependent diffusivity
g = 2+sin(x);
gx = rx.*(Dr*g);

%set lifting function to go from homogenous to inhomogenous solution
theta = a*x + ((b-a)/2/(xR-xL))*(x-xL).^2;
thetax = a + ((b-a)/(xR-xL))*(x-xL);
thetaxx = (b-a)/(xR-xL);

rhs = -cos(x).*(2*sin(x)+3); %rhs of inhomog problem
stuff = rhs +theta - g.*thetaxx - gx.*thetax;
f = J.*((invV'*invV)*stuff);

[A] = HelmCstab1Dsdd(g);

% Solve (homogenous) Problem
solvec = A\f(:);
v = reshape(solvec,Np,K);

%recover imhomogenous solution
u = v + theta;

uexact = cos(x);

figure(2);
plot(x(:),u(:),x(:),uexact(:));

norm(u(:)-uexact(:),2)
