function hatGrads = getHatGrads(mesh)
%GETHATGRADS    Gradients of the hat functions.
%   HATGRADS = GETHATGRADS(MESH) returns the gradients of the hat functions
%   (nodal basis of the elementwise linear, globally continuous functions)
%   on the simplicial mesh MESH. HATGRADS is a (DIM+1)-cell array, so that
%   HATGRADS{j}(i,:) is the gradient of the hat function corresponding to
%   the node MESH.ELEMENTS(i,j) on the element MESH.ELEMENTS(i,:).
%
%   Works for arbitrary-dimensional meshes.
%
%   Author: Josef Kemetmueller - 16.12.2013

if isfield(mesh,'hatGrads')
    hatGrads = mesh.hatGrads;
    return;
end
assert(dimMesh(mesh)<=dimSpace(mesh), ...
                ['Your %dD-simplices are lying in %dD-space.' ...
                'I''m not sure there even is a ''gradient''.'], ...
                dimMesh(mesh), dimSpace(mesh));

[X,hatGrads] = deal(cell(1,dimMesh(mesh)+1));
for d = 1:dimMesh(mesh)+1
    X{d} = mesh.coordinates(mesh.elements(:,d),:);
end
% "Modulo operation" m(i)=^=mod(i,dimMesh(mesh)+1) for small i,
%  with 0 =^= dimMesh(mesh)+1.
mo = [1:dimMesh(mesh)+1, 1:dimMesh(mesh)+1];
if dimMesh(mesh)==1
    divideByLengthSquared = @(X) bsxfun(@rdivide,X,dot(X,X,2));
    for i = 1:dimMesh(mesh)+1
        hatGrads{i} = divideByLengthSquared(X{mo(i+1)}-X{mo(i)});
    end
elseif dimMesh(mesh)==2 && dimSpace(mesh)<=3
    switch dimSpace(mesh)
        % orth(Xa,Xb) generates a vector orthogonal to the edge
        % [Xa,Xb] lying inside the triangle.
        case 2 %2D mesh embedded in R^2
            orth = @(Xa,Xb) [-Xa(:,2)+Xb(:,2), Xa(:,1)-Xb(:,1)];
        case 3 %2D mesh embedded in R^3
            N_3D = cross(X{3}-X{1},X{2}-X{1},2);
            orth = @(Xa,Xb) cross(N_3D,Xa-Xb,2);
    end
    for i = 1:dimMesh(mesh)+1
        Ni = orth(X{mo(i+2)},X{mo(i+1)});
        hatGrads{i} = bsxfun(@rdivide,Ni,dot(Ni,X{i}-X{mo(i+1)},2));
    end
elseif dimMesh(mesh)==3 && dimSpace(mesh)==3
    % Gradient has the same direction as normal vector of opposite
    % face. The length is 1/(length point to opposite face).
    for i = 1:dimMesh(mesh)+1
        Ni = cross(X{mo(i+2)}-X{mo(i+1)}, X{mo(i+3)}-X{mo(i+1)},2);
        hatGrads{i} = bsxfun(@rdivide,Ni,dot(Ni,X{i}-X{mo(i+1)},2));
    end
% Feel free to add the 4D case explicitly.
% elseif dimMesh(mesh)==4 && dimSpace(mesh)==4
else % general case for dimMesh <= dimSpace
    nE = numElements(mesh);
    [hatGrads{:}] = deal(zeros(nE, dimSpace(mesh)));
    for el = 1:nE % Well, you have to pay top dollar for R^n...
        nodes = mesh.elements(el,:);
        % MATLABs least squares solution generates a valid gradient vector
        % for d<n too, using the following approach! (This wouldn't work if
        % we were to consider the extended system: [1,1,1;coord]\[0,0;1,0;0,1] !!!)
        B = bsxfun(@minus, mesh.coordinates(nodes(2:end), :), ...
                           mesh.coordinates(nodes(1), :) ).';
        grads = B\eye(dimSpace(mesh));
        hatGrads{1}(el,:) = -sum(grads,1);
        for i = 2:dimMesh(mesh)+1
            hatGrads{i}(el,:) = grads(i-1,:);
        end
    end
end
end
