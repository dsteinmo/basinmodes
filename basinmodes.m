load(MESH_FILE);


StartUp2D;  %had to tweak NODETOL in this - apparently 1e-12 is too small.
BuildBCMaps2D;

for ii=1:NUM_H_REFINES
    refineflag= ones(K,1);
    Hrefine2D(refineflag);
    StartUp2D;
end

%interpolate bathymetry profile to unstructured mesh.
if ANALYTIC_DEPTH == false
    H = interp2(depthdata.x,depthdata.y,depthdata.H,x,y,'linear');
else
    H = H_analytic(x, y);
end

if ~isempty(find(vmapP==0, 1))
    disp('vmapP is bad node map, try relaxing NODETOL in StartUp2D');
    return;
end

A = dgint(ones(Np,K),V,J); %first, get the area of our basin by integrating '1'.

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
[OpNeu,MM] = PoissonIPDG2Dsdd_eig(h);

[efuns,d,convflag]=eigs(OpNeu,MM,numpot+1,'SM');   

if convflag ~= 0
    disp('eigenvalues didn''t converged! (normal laplacian operator)');
end

d= diag(d);
[lambda,inds] = sort(d,'ascend');
efuns = efuns(:,inds);

%get rid of mode zero
lambda = lambda(2:end);
efuns = efuns(:,2:end);

phi = cell(numpot,1);  %Allocate cell array for homogenous Neuman eigenfns

disp('checking normalization of phis...');
for jj=1:numpot
    phitmp = efuns(:,jj);
    
    phi{jj} = reshape(phitmp,Np,K);

    %normalize eigenfunction
    alpha = dgint(phi{jj}.*phi{jj},V,J);
    phi{jj} = phi{jj}/sqrt(alpha);
    phi{jj} = phi{jj}*(sqrt(A)*c0*Hbar)/sqrt(lambda(jj));
    
    disp(['Rel. error in normalization: ' num2str((lambda(jj)*dgint(phi{jj}.*phi{jj},V,J)-A*c02*Hbar^2)/(A*c02*Hbar^2))]); %works good
    
end
disp('done check.');

%Build Laplacian operator with Dirichlet BC's 
ids = find(BCType == Wall | BCType == Neuman);
BCType(ids) = Dirichlet;  

[OpDir,MM] = PoissonIPDG2Dsdd_eig(hinv);  

[efuns,d,convflag] = eigs(OpDir,MM,numstrm,'SM');   
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
end

%Now construct irrotational and solenoidal basis functions
%vphi_j = - grad phi_j, and vpsi_j = hinv (khat cross grad psi_j) .
vphix = cell(numpot,1);
vphiy = cell(numpot,1);

dphi = zeros(3*Nfp,K);
for jj = 1:numpot
    phitmp = phi{jj};
    [phi_x,phi_y] = Grad2D(phitmp);
    
    %Surface integral contributions to gradient. 
    %(Not sure if this is necessarily needed or not)
    dphi(:) = phitmp(vmapM) - phitmp(vmapP);
    
    fluxphix =nx.*(dphi/2); fluxphiy =ny.*(dphi/2);
    
    vphix{jj} = (phi_x - LIFT*(Fscale.*fluxphix));
    vphiy{jj} = (phi_y - LIFT*(Fscale.*fluxphiy)); 
end

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

%Compute constant coefficients (as in Rao&Schwab 1976)
nu = sqrt(c02.*lambda);
a = zeros(numpot,numpot); b = zeros(numpot,numstrm);
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
disp('done.');


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
    PlotField2D(N,x,y,real(eta{jj+start})); 
    view([0 90]); 
    axis equal;
    axis tight;
     
    colorbar; 
    title(['T= ' num2str(2*pi/abs(myeigs(jj+start))/3600) ' h']); 
    drawnow;
end

if DUMP_TO_FILE == true
    save('basinmodes_output.mat');
end