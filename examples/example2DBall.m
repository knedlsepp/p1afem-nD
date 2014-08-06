%% Just some random data
ballmesh = ball(100);
sphereWithHole = cleanMesh(genMesh(ballmesh.bd(1).elements(2:end,:),ballmesh.coordinates,'any'));
mesh = genBisectionMesh(sphereWithHole);
mesh.bd(2).elements = []; % Set neumann boundary to empty.
f = @(X) -sum(abs(X-2),2).^(1/2);
%%
g = @(X) zeros(size(X,1),1); % No boundary!
uD = @(X) zeros(size(X,1),1);
nEmax = 1e5;
rho = 0.1;
discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD, 2);
quadDeg = 2;

[x,energy] = solvePoisson(mesh, f, g, uD, discretizer, quadDeg);
surfMesh(mesh, x);
[x, mesh, indicators] = adaptivePoisson(mesh, f, g, uD, nEmax, rho, ...
                                        discretizer, quadDeg, @surfSolution);