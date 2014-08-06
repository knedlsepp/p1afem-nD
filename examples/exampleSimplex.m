%% Just some random data
mesh = genBisectionMesh(genMesh(1:4, [zeros(1,3);eye(3)], 1:3));
f = @(X) -sum(abs(X),2).^(1/2);
g = @(X) 1*X(:,2).^3;
uD = @(X) 1/2*X(:,1);
nEmax = 1e5;
rho = 0.1;
discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD, 2);
quadDeg = 2;

[x, mesh, indicators] = adaptivePoisson(mesh, f, g, uD, nEmax, rho, ...
                                        discretizer, quadDeg, @surfSolution);