function S = inner_L2_gradP1_gradP1_L2(mesh, s, quadDeg)
%INNER_L2_GRADP1_GRADP1_L2 L2-inner product of a scalar function s and the
%gradients of two elementwise linear, globally continuous functions.
%
%   M = inner_L2_gradP1_gradP1_L2(MESH, s) returns the matrix that corresponds
%   to the variable coefficient Poisson equation:
%       \int s*grad(u)*grad(v) dx
%
%   Author: Josef Kemetmueller - 23.05.2015
if ~isfield(mesh,'volumes')
    mesh.volumes = getElementVolumes(mesh);
end
if ~exist('quadDeg','var');
    quadDeg = 3;
end
nC = numCoordinates(mesh);
nE = numElements(mesh);

% (dimMesh+1) hat-functions per Element => nEx(dimMesh+1)x(dimMesh+1) contributions.

hatGrads = getHatGrads(mesh);
[S,I,J] = deal(zeros(nE,dimMesh(mesh)+1,dimMesh(mesh)+1));

sP0 = L2ProjectL2ToP0(mesh, s, quadDeg);
for i = 1:dimMesh(mesh)+1
    for j = 1:dimMesh(mesh)+1
		% This works because elementwise P0_P0 scalar product can be interpreted as nE local P0 functions.
		% For higher order elements this is not that simple. See e.g. inner_P1_P1_L2
        % The variable coefficient can simply be projected to P0 space and
        % multiplied with the other elementwise-constant function.
        S(:,i,j) = sP0.*inner_P0_P0_L2(mesh,hatGrads{i},hatGrads{j},'elementwise');
        I(:,i,j) = mesh.elements(:,i);
        J(:,i,j) = mesh.elements(:,j);
    end
end
%% Compute the matrix
    S = sparse(I(:),J(:),S(:),nC,nC);
end
