%% You have to download p1afem2d to make this work
load elements.dat
load coordinates.dat
load dirichlet.dat
load neumann.dat

mesh = genBisectionMesh(genMesh(elements, coordinates, dirichlet, neumann), 'forceRefine');

f = @(X) ones(size(X,1),1);
qeq = @(a,b) abs(a-b)<1e-9*abs(a); % quite equal
g = @(X) double((qeq(X(:,1),-1) & X(:,2)>=0) | (qeq(X(:,2),1) & X(:,1)>=0));
uD = @(X) zeros(size(X,1),1);

nEmax = 1e4;
rho = 0.1;
discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD, 2);
quadDeg = 2;

[x, mesh, indicators] = adaptivePoisson(mesh, f, g, uD, nEmax, rho, ...
                                        discretizer, quadDeg, @surfSolution);