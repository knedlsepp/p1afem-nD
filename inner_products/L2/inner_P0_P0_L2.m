function S = inner_P0_P0_L2(mesh, P0_1, P0_2, dimPhi)
%INNER_P0_P0_L2    L2-inner product of two piecewise constant functions.
%   S = INNER_P0_P0_L2(MESH, P0_1, P0_2) returns the inner_product (P0_1,P0_2)_L2.
%
%   V = INNER_P0_P0_L2(MESH, P0_1) returns the vector V of P0_1 tested with
%   all basis functions V(i) == INNER_P0_P0_L2(MESH, P0_1, CHI(i)), where
%   CHI(i) is the basis function which is one on the element(i,:) and zero
%   everywhere else.
%
%   M = INNER_P0_P0_L2(MESH) or M = INNER_P0_P0_L2(MESH, [], [], DIMP0_1)
%   returns the matrix M(i,j) = INNER_P0_P0_L2(MESH, CHI(i), CHI(j)), where
%   CHI(i) is again the basis function which is one on the element(i,:) and
%   zero everywhere else.
%   [The so called 'mass matrix' for P^0 functions].
%   Author: Josef Kemetmueller - 16.12.2013
if ~isfield(mesh,'volumes')
    mesh.volumes = getElementVolumes(mesh);
end
if exist('dimPhi','var') && strcmpi(dimPhi, 'elementwise')
    %% Matrix reuse:
    % S = dot(P0_1, inner_P0_P0_L2_mat(mesh, 1)*P0_2, 2);
    %% More speed:
    S = mesh.volumes.*dot(P0_1, P0_2, 2);
elseif ~exist('P0_1','var')
    S = inner_P0_P0_L2_mat(mesh, 1);
elseif ~exist('P0_2','var')
    S = inner_P0_P0_L2(mesh, P0_1, []);
elseif isempty(P0_1) && isempty(P0_2)
    S = inner_P0_P0_L2_mat(mesh, dimPhi);
elseif isempty(P0_1)
    S = inner_P0_P0_L2_mat(mesh, 1)*P0_2;
elseif isempty(P0_2)
    S = inner_P0_P0_L2_mat(mesh, 1)*P0_1;
else
    S = sum(dot(P0_1, inner_P0_P0_L2_mat(mesh, 1)*P0_2, 2));
end
end

function S = inner_P0_P0_L2_mat(mesh, dim)
nE = numElements(mesh);
S = spdiags(repmat(mesh.volumes(:),dim,1),0,dim*nE,dim*nE);
end
