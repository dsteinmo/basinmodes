function uapprox = mypreconDG2D(RHSfine,L,U,P,Q,InterpMat,CoarseMat,Npfine,Npcoarse,K)

    f = reshape(RHSfine,Npfine,K);
    
    %coarsen the RHS
    fcoarse = CoarseMat*f;
    
    %Solve the coarse problem
    ucoarse = Q*(U\(L\(P*fcoarse(:))));
    
    ucoarse = reshape(ucoarse,Npcoarse,K);
    
    %Smooth coarse solution and take it as
    %an approximation to the fine solution
    uapprox = InterpMat*ucoarse;
    
    %return vector form
    uapprox = uapprox(:);
end