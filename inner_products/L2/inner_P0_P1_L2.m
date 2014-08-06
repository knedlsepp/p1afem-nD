function S = inner_P0_P1_L2(mesh, P0, P1, dimPhi)
%INNER_P0_P1_L2    L2-inner product of an elementwise constant and an
%elementwise linear, globally continuous function.
%   S = INNER_P0_P1_L2(MESH, P0, P1) returns the inner product (P0,P1)_L2.
%
%   V = INNER_P0_P1_L2(MESH, P0) returns the vector V of P0 tested with
%   all hat functions V(i) == INNER_P0_P1_L2(MESH, P0, HAT(i))
%
%   V = INNER_P0_P1_L2(MESH, [], P1) returns the vector V of P1 tested with
%   all basis functions CHI(i): V(i) == INNER_P0_P1_L2(MESH, CHI(i), P1),
%   where CHI(i) is the basis function which is one on the element(i,:) and
%   zero everywhere else.
%
%   M = INNER_P0_P1_L2(MESH) or M = INNER_P0_P1_L2(MESH, [], [], DIMPHI)
%   returns the matrix M(i,j) = INNER_P0_P1_L2(MESH, CHI(i), HAT(j)).
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~exist('P0','var')
    S = inner_P0_P1_L2_mat(mesh, 1);
elseif ~exist('P1','var')
    S = inner_P0_P1_L2(mesh, P0, []);
elseif isempty(P0) && isempty(P1)
    S = inner_P0_P1_L2_mat(mesh, dimPhi);
elseif isempty(P0)
    S = inner_P0_P1_L2_mat(mesh, 1)*P1;
elseif isempty(P1)
    S = inner_P0_P1_L2_mat(mesh, 1)'*P0;
elseif exist('dimPhi','var') && strcmpi(dimPhi, 'elementwise')
    S = dot(P0, inner_P0_P1_L2_mat(mesh, 1)*P1, 2);
else
    S = sum(dot(P0, inner_P0_P1_L2_mat(mesh, 1)*P1, 2));
end
end

function S = inner_P0_P1_L2_mat(mesh, dim)

nC = numCoordinates(mesh);
nE = numElements(mesh);

S = repmat(mesh.volumes/(dimMesh(mesh)+1), 1, dimMesh(mesh)+1);
I = repmat(1:nE, 1, dimMesh(mesh)+1);
J = mesh.elements;

if dim==1
    S = sparse(I(:),J(:),S(:),nE,nC);
else
    S = sparseBlockDiag(I,J,S,nE,nC,dim);
end


end


