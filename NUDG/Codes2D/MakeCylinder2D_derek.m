function MakeCylinder2D_derek(faces, ra,xo,yo)

% Function: MakeCylinder2D(faces, ra, xo, yo)
% Purpose:  Use Gordon-Hall blending with an isoparametric map to modify a list
%           of faces so they conform to a cylinder of radius r centered at (xo,yo)
Globals2D;

NCurveFaces = size(faces,1);
vflag = zeros(size(VX));
figure(1); clf; PlotMesh2D; hold on;
for n=1:NCurveFaces 

  % move vertices of faces to be curved onto circle
  k = faces(n,1)
  f = faces(n,2)
  v1 = EToV(k, f); v2 = EToV(k, mod(f,Nfaces)+1);

  % compute polar angles of start and end face vertices relative to circle center
  theta1 = atan2(VY(v1)-yo,VX(v1)-xo); 
  theta2 = atan2(VY(v2)-yo,VX(v2)-xo);

  % move vertices onto circle
  newx1 = xo + ra*cos(theta1); newy1 = yo + ra*sin(theta1);
  newx2 = xo + ra*cos(theta2); newy2 = yo + ra*sin(theta2);

  disp('old vertices');
  disp(['v1 (' num2str(VX(v1)) ',' num2str(VY(v1)) ')']);
  disp(['v2 (' num2str(VX(v2)) ',' num2str(VY(v2)) ')']);
  
  plot(VX(v1),VY(v1),'.g');
  plot(VX(v2),VY(v2),'.g');
  
  % update mesh vertex locations
  VX(v1) = newx1; VX(v2) = newx2; VY(v1) = newy1; VY(v2) = newy2; 
  
  disp('new vertices');
  disp(['v1 (' num2str(VX(v1)) ',' num2str(VY(v1)) ')']);
  disp(['v2 (' num2str(VX(v2)) ',' num2str(VY(v2)) ')']);
  
  plot(VX(v1),VY(v1),'*r');
  plot(VX(v2),VY(v2),'*r');

  % store modified vertex numbers
  vflag(v1) = 1;  vflag(v2) = 1;
  
  %pause;
end

%keyboard;

% map modified vertex flag to each element
vflag = vflag(EToV);

% locate elements with at least one modified vertex
ks = find(sum(vflag,2)>0);

% build coordinates of all the corrected nodes
va = EToV(ks,1)'; vb = EToV(ks,2)'; vc = EToV(ks,3)';
x(:,ks) = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y(:,ks) = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

plot(x(:,ks),y(:,ks),'b.')
%this is still straight-sided, just making sure that
%triangle vertices are shifted onto the curvilinear boundary, and adjust
%nodes for these triangles accordingly

for n=1:NCurveFaces  % deform specified faces
  k = faces(n,1); f = faces(n,2);

  % find vertex locations for this face and tangential coordinate
  if(f==1) v1 = EToV(k,1); v2 = EToV(k,2); vr = r; end
  if(f==2) v1 = EToV(k,2); v2 = EToV(k,3); vr = s; end
  if(f==3) v1 = EToV(k,1); v2 = EToV(k,3); vr = s; end
  fr = vr(Fmask(:,f));
  x1 = VX(v1); y1 = VY(v1); x2 = VX(v2); y2 = VY(v2);

  %Derek: get end-points of the curvilinear face, then
  %parameterize it.
  
  %first get end-points of parameterization (theta1, theta2)
  theta1 = atan2(y1-yo, x1-xo); theta2 = atan2(y2-yo, x2-xo);

  % check to make sure they are in the same quadrant
  if ((theta2 > 0) & (theta1 < 0)), theta1 = theta1 + 2*pi; end;
  if ((theta1 > 0) & (theta2 < 0)), theta2 = theta2 + 2*pi; end;
  
  % Distribute N+1 face nodes by arc-length along edge,
  
  % Derek: this basically gives your parameter (theta) an LGL
  %spacing between the two end-points [theta1,theta2].
  %the formula should generalize for arbitrary arclength parameterizations.
  theta = 0.5*theta1*(1-fr) + 0.5*theta2*(1+fr);

  %keyboard;
  % evaluate deformation of coordinates (along face)
  % Derek: basically (xnew) - xold
  % where xnew is evaluated using the parameterization
  fdx = (xo + ra*cos(theta))-x(Fmask(:,f),k); 
  fdy = (yo + ra*sin(theta))-y(Fmask(:,f),k);
  
  % build 1D Vandermonde matrix for face nodes and volume nodes
  Vface = Vandermonde1D(N, fr);  Vvol  = Vandermonde1D(N, vr);
  % compute unblended volume deformations 
  vdx = Vvol*(Vface\fdx); vdy = Vvol*(Vface\fdy);

  %blending stuff should all stay the same.
  
  % blend deformation and increment node coordinates
  ids = find(abs(1-vr)>1e-7); % warp and blend
  if(f==1) blend = -(r(ids)+s(ids))./(1-vr(ids)); end;
  if(f==2) blend =      +(r(ids)+1)./(1-vr(ids)); end;
  if(f==3) blend = -(r(ids)+s(ids))./(1-vr(ids)); end;

  figure(31);
  plot(x(:,k),y(:,k),'.-');
  title('pre-blend');
  drawnow;
  
  %this is where the actual "curving" takes place.
  x(ids,k) = x(ids,k)+blend.*vdx(ids);
  y(ids,k) = y(ids,k)+blend.*vdy(ids);
  
  figure(32);
  plot(x(:,k),y(:,k),'.-');
  title('post-blend');
  drawnow;
  
  pause;
end

% repair other coordinate dependent information
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);
[rx,sx,ry,sy,J] = GeometricFactors2D(x, y,Dr,Ds);
[nx, ny, sJ] = Normals2D(); Fscale = sJ./(J(Fmask,:));
return
