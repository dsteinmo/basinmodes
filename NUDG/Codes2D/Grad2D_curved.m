function [ux,uy] = Grad2D_curved(u)
% function [ux,uy] = Grad2D(u);
% Purpose: Compute 2D gradient field of scalar u
%
% Uses nodal evaluation on straight triangles and cubature 
% integration for de-aliasing on curved triangles.

    global Dr Ds rx sx ry sy straight curved Np K cub;

    %allocate needed arrays
    ux = zeros(Np,K); uy = ux;

    %straight-sided elements first (nodal):
    ur_str = Dr*u(:,straight);
    us_str = Ds*u(:,straight);

    ux(:,straight) = rx(:,straight).*ur_str + sx(:,straight).*us_str;
    uy(:,straight) = ry(:,straight).*ur_str + sy(:,straight).*us_str;

    %now deal with curved elements
    ur_cur = cub.Dr*u(:,curved);
    us_cur = cub.Ds*u(:,curved);

    ux(:,curved) = cub.V'*(cub.W(:,curved).*(  cub.sx(:,curved).*us_cur ...
                                             + cub.rx(:,curved).*ur_cur));

    uy(:,curved) = cub.V'*(cub.W(:,curved).*(  cub.sy(:,curved).*us_cur ...
                                             + cub.ry(:,curved).*ur_cur));

    %catch is that on curved elements we actually computed Stiffness-matrix
    %times u, need to multiply by specialized inverse mass-matrix.
    for m=1:length(curved)
        k = curved(m);
        mmCHOL = cub.mmCHOL(:,:,k);
        ux(:,k) = mmCHOL\(mmCHOL'\ux(:,k));
        uy(:,k) = mmCHOL\(mmCHOL'\uy(:,k));
    end

end