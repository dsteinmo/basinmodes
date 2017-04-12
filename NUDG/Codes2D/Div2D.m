function [divu] = Div2D(u,v)

% function [divu] = Div2D(u,v);
% Purpose: Compute the 2D divergence of the vector field (u,v)

%Globals2D;
global Dr Ds rx sx ry sy

ur = Dr*u; us = Ds*u; vr = Dr*v; vs = Ds*v;
divu = rx.*ur + sx.*us + ry.*vr + sy.*vs;
return;
