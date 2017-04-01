function curvedElements = MakeCurvedEdges_derek(BCFlag,xs,ys,t)

% Function: MakeCylinder2D(BCFlag, xs, ys, t)
% Purpose:  Use Gordon-Hall blending with an isoparametric map to modify a list
%           of faces so they conform to a spline boundary. Parameterized by
%           xs(t),ys(t)

%note: I think this will break for 'exant' corners.
%ie, if an element has two edges on a boundary. 
%should perhaps put in code to do nothing if we detect one of those
%if it starts to be a problem

Globals2D;

[k,f] = find(BCType==BCFlag);


%If an element has more than 1 face on the boundary, exclude
%one of the faces. Since keeping it breaks things (negative jacobian).
[k,uniqind] = unique(k,'first');
f = f(uniqind);

faces =[k,f];

%get a fine evaluation of our parametric spline, for pushing boundary nodes
%onto.

tt = linspace(min(t),max(t),100*length(t)); %100 fairly arbitrary, may not need to be that fine.
xsval = ppval(xs,tt);  %evaluate spline at
ysval = ppval(ys,tt);  %finer parameter space.

NCurveFaces = size(faces,1);
vflag = zeros(size(VX));

%diagnostic plot
%figure(1); clf; PlotMesh2D; hold on;

%stuff for skipping faces that don't appear to having corresponding spline
%points
skip=[];
%get domain length-scale
L=max(abs(x(:)));
W=max(abs(y(:)));
L=max([L W]);
mytol = 0.3*L;
for n=1:NCurveFaces 

  % move vertices of faces to be curved onto spline
  k = faces(n,1);
  f = faces(n,2);
  v1 = EToV(k, f); 
  v2 = EToV(k, mod(f,Nfaces)+1);

  %find minimum square distance from existing node to point on spline
  v1_dists2 = (VX(v1)-xsval).^2 + (VY(v1)-ysval).^2;
  v2_dists2 = (VX(v2)-xsval).^2 + (VY(v2)-ysval).^2;
  
  [dist1 v1s_ind] = min(v1_dists2);
  [dist2 v2s_ind] = min(v2_dists2);
  
  
  if dist1 > mytol || dist2 > mytol
      %we don't appear to have a spline to deform onto
      %so we'll do nothing here
      skip=[skip;n];
  else
      %if we found a spline to deform onto, then move triangle 
      %vertices around...
      
      %set new vertex coordinates to closest spline point coordinates
      newx1 = xsval(v1s_ind); newy1 = ysval(v1s_ind);
      newx2 = xsval(v2s_ind); newy2 = ysval(v2s_ind);

      %diagnostic plot
      %plot([VX(v1) VX(v2)],[VY(v1) VY(v2)],'.r');
      %plot([newx1 newx2],[newy1 newy2],'^g');

      %update mesh vertex locations
      VX(v1) = newx1; VX(v2) = newx2; VY(v1) = newy1; VY(v2) = newy2; 

      % store modified vertex numbers
      vflag(v1) = 1;  vflag(v2) = 1;
  end
  
  %pause;
end


% map modified vertex flag to each element
vflag = vflag(EToV);

% locate elements with at least one modified vertex
ks = find(sum(vflag,2)>0);

% build coordinates of all the corrected nodes
va = EToV(ks,1)'; vb = EToV(ks,2)'; vc = EToV(ks,3)';
x(:,ks) = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y(:,ks) = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

%plot(x(:,ks),y(:,ks),'b.')
%this is still straight-sided, just making sure that
%triangle vertices are shifted onto the curvilinear boundary, and adjust
%nodes for these triangles accordingly

%for n=1:NCurveFaces  % deform specified faces

curvedFaces = setdiff(1:NCurveFaces,skip);
for n= curvedFaces
  k = faces(n,1); f = faces(n,2);

  % find vertex locations for this face and tangential coordinate
  if(f==1) v1 = EToV(k,1); v2 = EToV(k,2); vr = r; end
  if(f==2) v1 = EToV(k,2); v2 = EToV(k,3); vr = s; end
  if(f==3) v1 = EToV(k,1); v2 = EToV(k,3); vr = s; end
  fr = vr(Fmask(:,f));
  x1 = VX(v1); y1 = VY(v1); x2 = VX(v2); y2 = VY(v2);

  %Derek: get end-points of the curvilinear face, then
  %parameterize it.
  
  %find minimum square distance from existing node to point on spline
  v1_dists2 = (x1-xsval).^2 + (y1-ysval).^2;
  v2_dists2 = (x2-xsval).^2 + (y2-ysval).^2;
  
  tol=NODETOL;
  v1s_inds = find(sqrt(v1_dists2) < tol);
  v2s_inds = find(sqrt(v2_dists2) < tol);
  
  %[dump v1s_ind] = min(v1_dists2);
  %[dump v2s_ind] = min(v2_dists2);
  
  if length(v1s_inds) > 1 && length(v2s_inds) > 1
      disp('two non-unique points on parameterized curved found, this shouldn''t happen.');
      error('bad parameterization.'); 
  elseif length(v1s_inds) == 1 && length(v2s_inds) == 1
      %the best case
      t1 = tt(v1s_inds);  %set end-points of parameter-space interval.
      t2 = tt(v2s_inds);
  elseif length(v1s_inds) == 1 && length(v2s_inds) == 2
      %disp('2 V2s match');
      %plot(xsval,ysval,'.-'); hold on;
      %plot(xsval(v2s_inds),ysval(v2s_inds),'*g');
      t1 = tt(v1s_inds);
      t2 = tt(v2s_inds);      
      lengths = abs(t1-t2);
      [dump correct_ind] = min(lengths);
      t2 = tt(v2s_inds(correct_ind));
  elseif length(v1s_inds) == 2 && length(v2s_inds) == 1
      %disp('2 V1s match');
      %plot(xsval,ysval,'.-'); hold on;
      %plot(xsval(v1s_inds),ysval(v1s_inds),'*r');
      %so choose the one with the smallest arclength      
      t1 = tt(v1s_inds);
      t2 = tt(v2s_inds);
      lengths = abs(t1-t2);
      [dump correct_ind] = min(lengths);
      t1 = tt(v1s_inds(correct_ind));
  %elseif catch other bad cases such as no matches?
  %should probably catch that in first loop above...
  else 
      disp('Some uncaught error has occurred due to non-unique points on parameterized curve.');
      error('bad parameterization.');
  end
      
  
  % Distribute N+1 face nodes by arc-length along edge,

  % Derek: this basically gives your parameter (t) an LGL
  %spacing between the two end-points [t1,t2].
  tLGL = 0.5*t1*(1-fr) + 0.5*t2*(1+fr);

  %equivalently (like we do with cheb), could write this as:
  %tLGL = ((fr+1)/2)*(t2-t1) + t1
  
  %keyboard;
  % evaluate deformation of coordinates (along face)
  % Derek: basically (xnew) - xold
  % where xnew is evaluated using the parameterization
  
  fdx = ppval(xs,tLGL) - x(Fmask(:,f),k);
  fdy = ppval(ys,tLGL) - y(Fmask(:,f),k);
  
  % build 1D Vandermonde matrix for face nodes and volume nodes
  Vface = Vandermonde1D(N, fr);  Vvol  = Vandermonde1D(N, vr);
  % compute unblended volume deformations 
  vdx = Vvol*(Vface\fdx); vdy = Vvol*(Vface\fdy);

  %blending stuff should all stay the same.
  %keyboard;
  
  
  % blend deformation and increment node coordinates
  ids = find(abs(1-vr)>1e-7); % warp and blend
  if(f==1) blend = -(r(ids)+s(ids))./(1-vr(ids)); end; %vr =r
  if(f==2) blend =      +(r(ids)+1)./(1-vr(ids)); end; %vr =s
  if(f==3) blend = -(r(ids)+s(ids))./(1-vr(ids)); end; %vr =s

  %myblend = zeros(Np,1);
  %myblend(ids) = blend;
  %figure(30);
  %subplot(3,1,1);
  %PlotField2D_1tri(N,r,s,myblend); view([0 90]); colorbar; drawnow;
  %subplot(3,1,2);
  %PlotField2D_1tri(N,r,s,vdx); view([0 90]); colorbar; drawnow;
  %subplot(3,1,3);
  %PlotField2D_1tri(N,r,s,vdy); view([0 90]); colorbar; drawnow;
  %figure(31);
  %PlotMesh2D;
  %hold on;
  %plot(x(:,k),y(:,k),'.-r');
  %title('pre-blend');
  %drawnow;
  %keyboard;
  
  %this is where the actual "curving" takes place.
  x(ids,k) = x(ids,k)+blend.*vdx(ids);
  y(ids,k) = y(ids,k)+blend.*vdy(ids);
  %Derek: will need to store the 'blend.*vdx(inds) factor and the ids 
  %for each curvilinear element to be able to map velocities to the
  %standard element for local post-processing.
  
  %plot(x(:,k),y(:,k),'^g');
  %title('post-blend');
  %drawnow;
  %hold off;
  %pause;
end

% repair other coordinate dependent information
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);
[rx,sx,ry,sy,J] = GeometricFactors2D(x, y,Dr,Ds);
[nx, ny, sJ] = Normals2D(); Fscale = sJ./(J(Fmask,:));

skiplen=length(skip);
disp(['Deformed ' num2str(NCurveFaces-skiplen) ' faces of ' num2str(NCurveFaces) ' possible boundary faces.']);

curvedElements = unique(faces(curvedFaces,1));

return
