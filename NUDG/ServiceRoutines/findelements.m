figure(7);

PlotMesh2D;
while true;
    disp('Zoom to region. Press any key when done.');
    pause;
    
    disp('select a vertex');
    [myx,myy] = ginput(1);
    mytol=50; %10m tolerance
    dists = sqrt((VX-myx).^2 + (VY-myy).^2);
    myV = find(dists<mytol);
    if length(myV) ~= 1
        disp('No vertex selected, assuming you''re done refining...');
        refineflag=0;
    else
        disp(['You clicked on vertex number: ' num2str(myV)]);
        myelement = mod(find(EToV==myV),K);
        if myelement == 0
            myelement = K;
        end
        
        disp('elements belonging to that vertex:');
        myelement
    end
        
end