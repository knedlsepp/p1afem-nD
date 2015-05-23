%% Solve the variable coefficient Poisson problem with point source RHS
%   div(c*grad(u)) = f+delta(x-x0)

%% Mesh
[X,Y] = ndgrid(linspace(-1,1,101));
P = [X(:),Y(:)];
DT = delaunayTriangulation(P);
mesh = genMesh(DT.ConnectivityList, P);
mesh.bd(1).elements = getBoundary(mesh); % Dirichlet boundary
mesh.bd(2).elements = []; % Neumann boundary
%% Setup
% x0 = [0,0]
pointSources = [0,0];
pointSourceValues = 1;
% c(x,y) = 1+sin(x)
c = @(X) 1+sin(X(:,1));
% f(x,y) = 2
f = @(X) 2*ones(size(X,1),1);
% Dirichlet boundary condition:
% uD(x,y) = 0
uD = @(X) zeros(size(X,1),1);
% Neumann boundary condition
% g(x,y) = 0
% (Neumann boundary: mesh.bd(2) chosen empty, so the value we choose won't matter)
g = @(X) zeros(size(X,1),1); % g = 0 

%%
discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD, 2);
quadDeg = 2;

%% 
dirichletBd = meshBd(mesh,1);
neumannBd = meshBd(mesh,2);
A = inner_L2_gradP1_gradP1_L2(mesh, c);
%% Compute right hand side.
b = inner_L2_P1_L2(mesh,      f, [], quadDeg) ...
  + inner_L2_P1_L2(neumannBd, g, [], quadDeg) ...
  + evalPointSources(mesh, pointSources, pointSourceValues) ...
  - A*discretizer(dirichletBd, uD);
%% Compute the discrete solution
x = discretizer(dirichletBd, uD);
freenodes = setdiff(1:numCoordinates(mesh), dirichletBd.elements(:));
x(freenodes,:) = A(freenodes,freenodes)\b(freenodes,:);

trisurf(mesh.elements,mesh.coordinates(:,1), mesh.coordinates(:,2), x,'facecolor','interp','edgecolor', 'none');
