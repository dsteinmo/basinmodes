[Nv,VX, VY, K, EToV,BCType,node,edge] = readmsh('ontario.msh');

save('ontario_gmsh.mat','Nv','VX','VY','K','EToV','BCType','node','edge','depthdata');
