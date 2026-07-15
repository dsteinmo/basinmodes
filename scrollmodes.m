%scrollmodes - cross-platform (MATLAB + GNU Octave)

jj = numstrm/2 + 1;
plotting = true;

hFig = figure('KeyPressFcn', @scrollmodes_keypress);

while plotting
    clf(hFig);
    pf2d(N,x,y,real(eta{jj*2-1})); colorbar;
    title(['modenum=' num2str(jj*2-1) '  \sigma/f=' num2str(scaledeigs(jj*2-1))]);

    setappdata(hFig,'lastkey','');   % reset before waiting
    uiwait(hFig);                    % blocks until uiresume() is called
    key = getappdata(hFig,'lastkey');

    if strcmp(key,'leftarrow') && (2*jj-1) ~= 1
        jj = jj-1;
    elseif strcmp(key,'rightarrow') && (2*jj-1) ~= nummodes
        jj = jj+1;
    else
        plotting = false;
    end
end

if ishandle(hFig)
    close(hFig);
end

