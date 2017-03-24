%'circular lake'modes data calculated using 200
%potentials and streamfunctions, with order 1 DG method on 2840 elements

%load data, if you haven't already: 
%load('circle_greatlake_modes.mat'); 

%Bottom is flat, so first 200 modes (the rotational modes, i.e., Rossby
%waves or shelf waves) are zero (to numerical precision). Gravity modes
%start at 201.

%mode=201; %gravest gravity/Kelvin mode = 201.
%mode=204;

%I believe half of the Kelvin modes (the negative frequency ones that are rotating clock-wise) are
%degenerate, and show up here as an artifice of the fact that the equations
%were set up such that the operator would be self-adjont (Hermitian).

%Note: mode 229 is the first super-inertial or Poincare mode. Why that is, I have no
%idea... i.e., What "portion" of the spectrum is allotted to the Kelvin modes?

%31 looks pretty good

figure(2); clf; colormap(darkjet);

%get angular frequency
freq = myeigs(mode);

T = 2*pi/abs(freq) %get period, so we'll plot one period.
numframes = 50;
t=0:T/numframes:T;  %make array of times to evaluate mode at.


for jj=1:length(t)
    pf2d(N,x,y,real(exp(1i*freq*t(jj))*eta{mode})); colorbar;
    %pf2d(N,x,y,real(exp(1i*rotmode_freq*t(jj))*rotmodeeta)); %colorbar;
    %caxis([-3 3])
    %caxis([-50 50]);
    view([0 90]); 
    axis tight; axis equal; colorbar; 
    %caxis([-4 4]);
    %title(['\sigma/f='  num2str(scaledeigs(mode))]);
    %title(['\sigma/f='  num2str(scaledeigs(mode)) ', t=' num2str(t(jj)) ' of ' num2str(T)]);
    title(['w/f=' sprintf('%1.2f, t=%1.2f of %1.2f',scaledeigs(mode),t(jj),T)]);
    print('-dpng',sprintf('frame%07d.png',jj));
    %t(jj)
    drawnow;
    %pause;
end
