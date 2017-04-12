function [u,time] = Heat1D(u,FinalTime)

% function [u] = Heat1D(u,FinalTime)
% Purpose  : Integrate 1D heat equation until 
%            FinalTime starting with initial condition, u.

Globals1D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np, K); 

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.25;dt   = CFL*(xmin)^2;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

%do this for LDG fluxes.
dt = 0.5*dt;

figure(1);

% outer time step loop 
for tstep=1:Nsteps
  for INTRK = 1:5
    timelocal = time + rk4c(INTRK)*dt;        

    % compute right hand side of 1D advection equations
    %[rhsu] = HeatCRHS1D(u,timelocal); %Central flux
    [rhsu] = HeatLDGRHS1D(u,timelocal); %LDG flux

    % initiate and increment Runge-Kutta residuals
    resu = rk4a(INTRK)*resu + dt*rhsu;  
    
    % update fields
    u = u+rk4b(INTRK)*resu;
    
    plot(x(:),u(:));
    title(['t=' num2str(time)]);
    ylim([-1 1]);
    drawnow;
  end;
  % Increment time
  time = time+dt;
end
return
