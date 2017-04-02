MSH_FILE = 'ontario2.msh';
MAT_FILE = 'ontario_gmsh2.mat';

[Nv,VX, VY, K, EToV,BCType,node,edge] = readmsh(MSH_FILE);

save(MAT_FILE,'Nv','VX','VY','K','EToV','BCType','node','edge','depthdata');
