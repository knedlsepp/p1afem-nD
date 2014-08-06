function P0 = L2ProjectP1ToP0(mesh, P1)
%L2PROJECTP1TOP0 L2-Projection from P1 to P0.
%   P0 = L2PROJECTP1TOP0(MESH, P1) returns the L2-projection onto the space
%   of elementwise constant functions.
%
%   Author: Josef Kemetmueller - 16.12.2013
nE = numElements(mesh);

% As this projection is just the elementwise mean, we compute it explicitly 
P0 = sum(reshape(P1(mesh.elements,:)',[],nE,dimMesh(mesh)+1),3)'/(dimMesh(mesh)+1);

% This would be an alternative to compute it, but since its a diagonal
% matrix, we do it explicitly:
%P0 = inner_P0_P0_L2(mesh)\inner_P0_P1_L2(mesh, [], P1);
