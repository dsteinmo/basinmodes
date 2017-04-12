function divu = Div2D_curved(u,v)
    global Dr Ds rx sx ry sy straight curved Np K cub;
   
    divu = zeros(Np,K);
    
    ur_str = Dr*u(:,straight); 
    us_str = Ds*u(:,straight); 
    vr_str = Dr*v(:,straight);
    vs_str = Ds*v(:,straight); 
    
    divu(:,straight) = rx(:,straight).*ur_str + sx(:,straight).*us_str ...
                     + ry(:,straight).*vr_str + sy(:,straight).*vs_str;
    
    %now deal with curved elements
    ur_cur = cub.Dr*u(:,curved); 
    us_cur = cub.Ds*u(:,curved); 
    vr_cur = cub.Dr*v(:,curved);
    vs_cur = cub.Ds*v(:,curved);
    
    %put in ux contribution
    %keyboard;
    divu(:,curved) = cub.V'*(cub.W(:,curved).*(  cub.sx(:,curved).*us_cur ...
                                             + cub.rx(:,curved).*ur_cur));
    %put in vy contribution                                         
    divu(:,curved) = divu(:,curved) + ...
                     cub.V'*(cub.W(:,curved).*(  cub.sy(:,curved).*vs_cur ...
                                             + cub.ry(:,curved).*vr_cur));
                                         
    %catch is that on curved elements we actually computed Stiffness-matrix
    %times u, need to multiply by specialized inverse mass-matrix.
    for m=1:length(curved)
        k = curved(m);
        mmCHOL = cub.mmCHOL(:,:,k);
        divu(:,k) = mmCHOL\(mmCHOL'\divu(:,k));
    end
   
end