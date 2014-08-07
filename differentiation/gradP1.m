function P0 = gradP1(mesh, P1)
%GRADP1    Gradient of elementwise linear, globally continuous function.
%   P0 = GRADP1(mesh, P1) returns the gradient of the function P1
%   represented via the nodal basis as a nC-by-dimP1 matrix, where nC is
%   the number of nodes in the mesh and dimP1 is the dimension of the
%   codomain of the function P1. P0 is a nE-by-dimSpace*dimP1 array
%   containing the gradients on all elements.
%
%   JacTr = reshape(P0(i,:), dimMesh, dimP1) is the transposed
%   jacobian-matrix on the i-th element.
%
%   Works for arbitrary-dimensional meshes.
%
%   Author: Josef Kemetmueller - 16.12.2013
nC = numCoordinates(mesh);
nE = numElements(mesh);
dimP1 = size(P1,2);
%
assert(size(P1,1)==nC, 'Function must be given as a P1 function');
%
P0 = zeros(nE,dimSpace(mesh),dimP1);
%%
hatGrads = getHatGrads(mesh);
for i = 1:dimP1
    for j = 1:dimMesh(mesh)+1
        % P0(T,:,i) is the gradient of the i-th component of P1 on the element T.
        P0(:,:,i) = P0(:,:,i) + bsxfun(@times, hatGrads{j}, P1(mesh.elements(:,j),i));
    end
end
P0 = reshape(P0,nE,[]);
end