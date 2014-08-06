%% Fichera cube with a point singularity at zero of the normal derivative.
load mesh_ficherasCube
mesh = genBisectionMesh(mesh);

f  = @(C) -(3/4)*(C(:,1).^2+C(:,2).^2+C(:,3).^2).^(-3/4);
uD = @(C) (C(:,1).^2+C(:,2).^2+C(:,3).^2).^(1/4);
g =  @(C) (1/2)*any(abs(C-1)<eps | abs(C+1)<eps,2).*(dot(C,C,2)).^(-3/4);
nEmax = 4e4;
rho = 0.1;
discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD, 2);
quadDeg = 2;

[x, mesh, indicators] = adaptivePoisson(mesh, f, g, uD, nEmax, rho, ...
                                        discretizer, quadDeg, ...
                                        @(mesh,x,I) surfSolution(mesh,x,I));