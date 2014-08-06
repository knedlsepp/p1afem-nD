%% Just some random data
mesh = genBisectionMesh(genMesh(1:2, [0,0;1,0], [1;2], []));

f = @(X) cos(X(:,1)*5*pi);
%%
g = @(X) -1e-3*ones(size(X,1),1);
uD = @(X) ones(size(X,1),1);
nEmax = 2e3;            
rho = 0.1;
discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD, 2);
quadDeg = 2;

[x,energy] = solvePoisson(mesh, f, g, uD, discretizer, quadDeg);
surfMesh(mesh, x);
[x, mesh, indicators] = adaptivePoisson(mesh, f, g, uD, nEmax, rho, ...
                                        discretizer, quadDeg, @surfSolution);