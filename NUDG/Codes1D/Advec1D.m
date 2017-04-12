function [u] = Advec1D(u, FinalTime)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;
time = 0;

%Let's try an integral, just for fun.
myfcn = sech(x-pi).^2;  %fcn we want to integrate
umodes = V\myfcn;  %transform nodes into modes
modezero = umodes(1,:);  %Get mode zero on each element
Po = 1/sqrt(2); %Scaling factor, since the height of order 0 ortho. legendre poly is 1/sqrt(2)
%modezero = V(1,:)*myfcn;
Klengths = x(Np,:) - x(1,:); %make an array of element lengths
myint = Po*(modezero*Klengths'); %dot product to sum up integral on each element

%exact = pi;
%exact = 2*exp(2*pi)*(-1 + exp(4*pi))./( (1+exp(2*pi))*(exp(2*pi) + exp(4*pi)) );

%error = myint - exact;


% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.75; dt   = CFL/(2*pi)*xmin; dt = .5*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

% advection speed
a = 2*pi;

% specifiy forcing function
force = 0*sech(2*(x-5)).^2;

% outer time step loop 
figure(1);
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u, timelocal, a);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
        u = u + dt*force; %add forcing term. could make time-dependent with an update step
        
        subplot(2,1,1);
        plot(x(:),u(:),'.-'); 
        axis([0 2*pi -1 1]);
        subplot(2,1,2);
        plot(x(vmapM),u(vmapM)-u(vmapP));
        %ylim([-0.1 1.2]);
        title(['t=' num2str(time)]);
        drawnow;
        
        
        %plot exact solution for sech^2 forcing?
        
    end
    % Increment time
    time = time+dt;
end
return
