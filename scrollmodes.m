%scrollmodes

jj=numstrm/2+1;
plotting = true;
while plotting==true
    pf2d(N,x,y,real(eta{jj*2-1})); colorbar;
    title(['modenum=' num2str(jj*2-1) '\sigma/f='  num2str(scaledeigs(jj*2-1))]);
    s=getkey('non-ascii');
    
    if strcmp(s,'leftarrow') && (2*jj-1) ~= 1
        jj=jj-1;
    elseif strcmp(s,'rightarrow') && (2*jj-1) ~= nummodes
        jj=jj+1;
    else
        plotting=false;
    end
end