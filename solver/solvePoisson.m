function [x, energy] = solvePoisson(mesh, f, g, uD, discretizer, quadDeg)
%SOLVEPOISSON    Solves the poisson problem with mixed Dirichlet-Neumann
%boundary conditions.
%   x = SOLVEPOISSON(MESH, f, g, uD, discretizer, quadDeg) solves the
%   Poisson problem:
%       - div(grad(u)) = f   in Omega
%                    u = uD  on the Dirichlet boundary
%               d/dn u = g   on the Neumann boundary,
%   using lowest order FEM, the equation resulting in
%
%      (grad u, grad v)_L2(Omega) = (f,v)_L2(Omega) + (g,v)_L2(Gamma_N),
%
%   with the side constraint u = uD on Gamma_D.
%   MESH is a simplicial mesh, generated using genMesh, where MESH.bd(1) is
%   the dirichlet boundary Gamma_D and MESH.bd(2) is the neumann boundary
%   Gamma_N, f/g/uD are all function handles, which can evaluate multiple
%   points at once, e.g. f(MESH.coordinates) must result in a
%   numCoordinates-by-1 array of the evaluations of
%   f(MESH.coordinates(j,:)). The values at the dirichlet boundary are
%   computed using a function handle DISCRETIZER of the signature
%   @(dirichletMesh, uD), e.g.
%       discretizer = @(dirichletBd, uD) L2ProjectL2ToP1(dirichletBd, uD);
%       discretizer = @(dirichletBd, uD) interpolateNodal(dirichletBd, uD);
%   The quadrature degree for the integration is given by quadDeg.
%
%   Works for arbitrary-dimensional meshes.
%
%   See also:
%   GENMESH
%   L2PROJECTL2TOP1
%   INTERPOLATENODAL
%   ADAPTIVEPOISSON
%
%   Author: Josef Kemetmueller - 28.03.2013

dirichletBd = meshBd(mesh,1);
neumannBd = meshBd(mesh,2);
nC = numCoordinates(mesh);
%% Compute stiffness matrix.
A = inner_gradP1_gradP1_L2(mesh);
%% Compute right hand side.
b = inner_L2_P1_L2(mesh,      f, [], quadDeg) ...
  + inner_L2_P1_L2(neumannBd, g, [], quadDeg) ...
  - A*discretizer(dirichletBd, uD);
%% Compute the discrete solution and the solution's energy norm.
x = discretizer(dirichletBd, uD);
freenodes = setdiff(1:nC, dirichletBd.elements(:));
x(freenodes,:) = A(freenodes,freenodes)\b(freenodes,:);
energy = sum(dot(x,A*x,1));