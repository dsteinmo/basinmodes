clear;
close all;

setPaths;

Globals2D;

% Polynomial order used for approximation 
N =1; 

%memory estimate: bytes=0; for i=1:length(blah), bytes = bytes + blah(i).bytes, end

%load mesh
load('ontario_gmsh.mat'); %disc with Csanady's radius choice = 67.5 km


StartUp2D;  %had to tweak NODETOL in this - apparently 1e-12 is too small.
BuildBCMaps2D;

%H-refinment, anyone?

numHRefines = 0;

for ii=1:numHRefines
    refineflag= ones(K,1);
    Hrefine2D(refineflag);
    StartUp2D;
end

USEMEANDEPTH = false; %use mean depth instead of full bathymetry?


numpot = 100;  %number of potential basis functions
numstrm = 100; %number of streamfunction basis functions


%interpolate bathymetry profile to unstructured mesh.
H = interp2(depthdata.x,depthdata.y,depthdata.H,x,y);

if ~isempty(find(vmapP==0, 1))
    disp('vmapP is bad node map, try relaxing NODETOL in StartUp2D');
    return;
end

%Physical params. (for model great lake, From Csanady 1967)
f= 2*(2*pi/(3600*24))*sin(43.7*pi/180); %43.7 deg. lattitude

g=9.81;


%Csanady's normal modes decomposition

% Hepi=12.5;        %Csanady calls "h",  epilimnion  thickness (m)
% Hhyp =27.5;       %Csanady calls "h'", hypolimnion thickness (m)
% 
% deltarho = 0.00201; %g/cm^3
% rhoprime = 1; %references density in g/cm^3
% 
% epsilon = deltarho/rhoprime; %relative density difference
% gprime = epsilon*g;  %reduced gravity
% 
% H1 = Hepi+Hhyp;  %total lake depth, (equivalent depth for barotropic mode)
% H2 = epsilon*((Hepi*Hhyp)/(Hepi+Hhyp)); %equivalent depth for internal mode.
%       %Note that Csanady absorb's "epsilon" into equivalent depth, while
%       %Kundu and others absorb epsilon into g, for a reduced gravity g'.
%       
% c1 = sqrt(g*H1);
% %c2 = sqrt(g*epsilon*((Hepi*Hhyp)/(Hepi+Hhyp)));
% c2 = sqrt(g*H2);

%H2 = epsilon


%g=9.81;
%H=1*ones(Np,K); %or set this to some function of (x,y)
%H = H2*ones(Np,K); %use Csanady's equivalent depth for baroclinic mode.

%Note: The final spectrum consists of (2n/3) gravitational modes
%and (n/3) rotational modes (modified by divergence). If H=const
%and f=const, the rotational modes vanish.

A = dgint(ones(Np,K),V,J); %first, get the area of our basin by integrating '1'.

%Define burger number Length scale = sqrt(area)
%S = c2/(f*sqrt(A));
%for circular basin, if S < (1/sqrt(2pi)) expect Kelvin modes, if S >
%(1/sqrt(2*pi)), then we don't.

Hbar = (1/A)*dgint(H,V,J); %basin mean depth

if USEMEANDEPTH == true
    H = Hbar*ones(Np,K);
end

c02 = g*Hbar;
c0 = sqrt(c02);

h = H/Hbar;  %non-dimensionalized depth
hinv = 1./h; %need this for the Dirichlet problem.


%Build Laplacian operator with Neuman BC's (assumes mesh has boundaries set
%to Wall or Neuman)
%[LapNeu,MM] = PoissonIPDG2D_eig(); %returns -lapacian (so pos. def.)
[OpNeu,MM] = PoissonIPDG2Dsdd_eig(h);

%%Second way of computing phis, get rid of mode zero.
[efuns,d,convflag]=eigs(OpNeu,MM,numpot+1,'SM');   
%[V,d]=eig(full(A));

if convflag ~= 0
    disp('eigenvalues didn''t converged! (normal laplacian operator)');
end

d= diag(d);
%oldd =d;
[lambda,inds] = sort(d,'ascend');
%[d,inds] = sort(d,'ascend');
efuns = efuns(:,inds);

%get rid of mode zero
lambda = lambda(2:end);
efuns = efuns(:,2:end);

phi = cell(numpot,1);  %Allocate cell array for homogenous Neuman eigenfns

disp('checking normalization of phis...');
for jj=1:numpot
    phitmp = efuns(:,jj);
    
    phi{jj} = reshape(phitmp,Np,K);
    
   % keyboard;
    
    %normalize eigenfunction
    alpha = dgint(phi{jj}.*phi{jj},V,J);
    phi{jj} = phi{jj}/sqrt(alpha);
    phi{jj} = phi{jj}*(sqrt(A)*c0*Hbar)/sqrt(lambda(jj));
    
    %phi{jj} = ((sqrt(A)*c0*Hbar)/sqrt(alpha)/sqrt(lambda(jj))).*phi{jj};
    
    %phi{jj} = ((sqrt(A)*c0*Hbar)/sqrt(alpha)/sqrt(lambda(jj))).*phi{jj};
    
    %what's the diff between next two lines?
    
    disp(['Rel. error in normalization: ' num2str((lambda(jj)*dgint(phi{jj}.*phi{jj},V,J)-A*c02*Hbar^2)/(A*c02*Hbar^2))]); %works good
    %disp(['Rel. error in normalization: ' num2str((dgint(lambda(jj)*phi{jj}.*phi{jj},V,J) -A*c02*Hbar^2)/(A*c02*Hbar^2))]); %works good
    
    
end
disp('done check.');
%%End second attempt at computing phis

%Build Laplacian operator with Dirichlet BC's 
ids = find(BCType == Wall | BCType == Neuman);
BCType(ids) = Dirichlet;  %change mesh's BC table entries to Dirichlet.

[OpDir,MM] = PoissonIPDG2Dsdd_eig(hinv);  %returns -lapacian (so pos. def.)

[efuns,d,convflag]=eigs(OpDir,MM,numstrm,'SM');   
if convflag ~= 0
    disp('eigenvalues didn''t converged! (normal laplacian operator)');
end

d= diag(d);
oldd =d;
[mu,inds] = sort(d,'ascend');
efuns = efuns(:,inds);


disp('checking normalization of psis...');
psi = cell(numstrm,1);  %Allocate cell array for homogenous Dirichlet eigenfns
for jj=1:numstrm
    psitmp = efuns(:,jj);
    psi{jj} = reshape(psitmp,Np,K); 
    
    %normalize eigenfunction
    alpha = dgint(psi{jj}.*psi{jj},V,J);
    psi{jj} = ((sqrt(A)*c0*Hbar)/sqrt(alpha)/sqrt(mu(jj))).*psi{jj};
    
    disp(['Rel. error in normalization: ' num2str((mu(jj)*dgint(psi{jj}.*psi{jj},V,J)-A*c02*Hbar^2)/(A*c02*Hbar^2))]); %works good
    %disp(['Relrao. error in normalization: ' num2str((mu(jj)*dgint(psi{jj}.*psi{jj},V,J)-A*c02*Hbar^2)/A*c02*Hbar^2)]);
end

%Now construct irrotational and solenoidal basis functions
%vphi_j = - grad phi_j, and
%vpsi_j = hinv (khat cross grad psi_j) .
vphix = cell(numpot,1);
vphiy = cell(numpot,1);

dphi = zeros(3*Nfp,K);
%disp('checking gradient of phi');
for jj = 1:numpot
    phitmp = phi{jj};
    %[vphix{jj},vphiy{jj}] = Grad2D(phitmp); %weak gradient
    [phi_x,phi_y] = Grad2D(phitmp);
    
    %Surface integral contributions to gradient. 
    %(Not sure if this is necessarily needed or not)
    dphi(:) = phitmp(vmapM) - phitmp(vmapP);
    
    fluxphix =nx.*(dphi/2); fluxphiy =ny.*(dphi/2);
    
    vphix{jj} = (phi_x - LIFT*(Fscale.*fluxphix));
    vphiy{jj} = (phi_y - LIFT*(Fscale.*fluxphiy)); 
end
%done irrotational

%now solenoidal:
%vpsix = hinv (-psi_y) , vpsiy = hinv (psi_x)
vpsix = cell(numstrm,1);
vpsiy = cell(numstrm,1);

dpsi = zeros(3*Nfp,K);
for jj=1:numstrm
    psitmp = psi{jj};
    [psi_x,psi_y] = Grad2D(psitmp);
    %Surface integral contributions to gradient. 
    dpsi(:) = psitmp(vmapM) - psitmp(vmapP);
    fluxpsix = nx.*(dpsi/2); fluxpsiy = ny.*(dpsi/2);
    
    psi_x = psi_x - LIFT*(Fscale.*fluxpsix);
    psi_y = psi_y - LIFT*(Fscale.*fluxpsiy);
    
    vpsix{jj} = psi_x;
    vpsiy{jj} = psi_y;
end
%done solenoidal basis functions.


%Compute constant coefficients (as in Rao&Schwab 1976)
nu = sqrt(c02.*lambda);
a = zeros(numpot,numpot); b = zeros(numpot,numstrm);
%c = zeros(numstrm,numpot); 
d = zeros(numstrm,numstrm);


%compute 'a' matrix:
for ii=1:numpot
    for jj=1:numpot
        aintegrand = -f.*(1/c02/A/(Hbar.^2)).*h.*(-vphix{ii}.*vphiy{jj} + vphiy{ii}.*vphix{jj});
        a(ii,jj) = dgint(aintegrand,V,J);
    end
end
%compute 'b' matrix:
for ii=1:numpot
    for jj=1:numstrm
        bintegrand = -f.*(1/c02/A/(Hbar.^2)).*(vphix{ii}.*vpsix{jj} + vphiy{ii}.*vpsiy{jj});
        b(ii,jj) = dgint(bintegrand,V,J);
    end
end
%get 'c' matrix, using the symmetry:
c = -b';
%compute 'd' matrix
for ii=1:numstrm
    for jj=1:numstrm
        dintegrand = f.*(1/c02/A/(Hbar.^2)).*hinv.*(-vpsiy{ii}.*vpsix{jj} + vpsiy{jj}.*vpsix{ii});
        d(ii,jj) = dgint(dintegrand,V,J);
    end
end

disp('checking coefficients symmetries');
norm(a+a',2)
norm(d+d',2)
norm(c+b',2)
norm(d+d',2)
disp('done that');


nudiag = diag(nu);

finalMat = 1i*[-a    -b   -nudiag;
                -c     -d   zeros(numstrm,numpot);
                nudiag  zeros(numpot,numstrm) zeros(numpot,numpot)];
       
%Sanity check: is finalMat hermitian? 
            
[efuns,myeigs] = eig(finalMat);
myeigs = diag(myeigs);
[dump,inds] = sort(abs(myeigs),'ascend');  %sort by absolute frequency.
myeigs = myeigs(inds);

periods = 2*pi./myeigs;
periods = periods/3600

scaledeigs = myeigs/f;

%Imma let chu finish, but...
%return;

efuns = efuns(:,inds);

nummodes=2*numpot+1*numstrm;

potential = cell(nummodes,1);
streamfcn = cell(nummodes,1);
eta = cell(nummodes,1);

for jj=1:nummodes
    potential{jj} = zeros(Np,K);
    eta{jj} = zeros(Np,K);
    for kk=1:numpot
        potential{jj} = potential{jj} + efuns(kk,jj)*phi{kk};
        eta{jj} = eta{jj} + efuns(kk+numpot+numstrm,jj)*(sqrt(lambda(kk))/c0)*phi{kk};
    end
    
end
for jj=1:nummodes
    streamfcn{jj} = zeros(Np,K);
    for kk=1:numstrm
        streamfcn{jj} = streamfcn{jj} + efuns(kk+numpot,jj)*psi{kk};
    end
end

%Plot the first 4 gravity modes.
figure(1); clf; colormap(darkjet);
start=numstrm+1;
for jj=1:2:8
    subplot(2,2,ceil(jj/2));
    PlotField2D(N,x,y,real(eta{jj+start})); view([0 90]); axis tight; colorbar; 
    %title(['\sigma/f='  num2str(abs(scaledeigs(jj+start)))]); 
    title(['T= ' num2str(2*pi/abs(myeigs(jj+start))/3600) ' h']); 
    caxis([-10 10]);
    drawnow;
end

save('ontario_allmodes.mat');
