%% Generate a mesh
mesh = genMesh([1,2,3;3,4,1], [0,0;1,0;1,1;0,1], [1,2], [2,3;3,4;4,1]);
mesh = genBisectionMesh(mesh);
%% Define poisson problem data
f = @(X) -sum(abs(X),2).^(1/2);
g = @(X) 1*X(:,2).^3;
uD = @(X) 1/2*X(:,1);
%% Define parameters for adaptive algorithm
nEmax = 1e5;
rho = 0.1;
discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD, 2);
quadDeg = 2;
%% Start computation
[x, refinedMesh, indicators] = adaptivePoisson(mesh, f, g, uD, nEmax, rho, ...
                                               discretizer, quadDeg, @surfSolution);
