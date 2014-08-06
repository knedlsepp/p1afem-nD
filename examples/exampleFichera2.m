%% Fichera cube with an edge singularity of the normal derivative.
load mesh_ficherasCube;
mesh = genBisectionMesh(mesh);

uD = @(C) ((C(:,1).^2 + C(:,2).^2 + C(:,3).^2).*((C(:,2).^2 + C(:,3).^2))).^(1/4);
f = @(C) -(C(:,1).^2+ 6*(C(:,2).^2+C(:,3).^2))./(4*((C(:,2).^2+C(:,3).^2).*(C(:,1).^2+C(:,2).^2+C(:,3).^2)).^(3/4));
g = @(C) (1/2)*any( abs(C(:,1)-1)<eps | abs(C(:,1)+1)<eps,2)...
                .*((C(:,2).^2+C(:,3).^2))./(((C(:,2).^2+C(:,3).^2).*(C(:,1).^2+C(:,2).^2+C(:,3).^2)).^(3/4)) ...
        + (1/2)*any(abs(C(:,2:3)-1)<eps | abs(C(:,2:3)+1)<eps,2)...
                .*((C(:,1).^2+2*(C(:,2).^2+C(:,3).^2)))./(((C(:,2).^2+C(:,3).^2).*(C(:,1).^2+C(:,2).^2+C(:,3).^2)).^(3/4));
nEmax = 1e5;
rho = 0.1;
discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD, 2);
quadDeg = 2;

[x, mesh, indicators] = adaptivePoisson(mesh, f, g, uD, nEmax, rho, ...
                                        discretizer, quadDeg, ...
                                        @(mesh,x,I) surfSolution(mesh,x,I));