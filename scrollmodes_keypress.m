function scrollmodes_keypress(src, event)
    setappdata(src,'lastkey',event.Key);
    uiresume(src);
end