function PlotDomain2D_black()

% function PlotDomain2D()
% Purpose: Show domain boundary

Globals2D;

hold on;
Nwall    = length(find(BCType(:)==Wall));
Ninflow  = length(find(BCType(:)==In));
Noutflow = length(find(BCType(:)==Out));
Ncyl     = length(find(BCType(:)==Cyl));


%ha = legend(types);
%set(ha, 'Fontsize', 16)

for k=1:K
  for f=1:Nfaces
    bc = BCType(k,f);
    ids = (k-1)*Nfp*Nfaces+(f-1)*Nfp+(1:Nfp);
    switch(bc)
      case Wall
	    plot3(Fx(ids), Fy(ids), 10+0*Fx(ids), 'k-');
      case Dirichlet
	   plot3(Fx(ids), Fy(ids), 10+0*Fx(ids), 'k-');
     case Neuman
       plot3(Fx(ids), Fy(ids), 10+0*Fx(ids), 'k-');  
      case In
	plot3(Fx(ids), Fy(ids), 10+0*Fx(ids), 'k--');
      case Out
	plot3(Fx(ids), Fy(ids), 10+0*Fx(ids), 'k:');
      case Cyl
	plot3(Fx(ids), Fy(ids), 10+0*Fx(ids),'k-.');
    end
  end
end

axis tight;
return;
