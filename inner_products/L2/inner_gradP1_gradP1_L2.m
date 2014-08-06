function S = inner_gradP1_gradP1_L2(mesh, phi, psi, dimPhi)
%INNER_GRADP1_GRADP1_L2 L2-inner product of the gradients of two
%elementwise linear, globally continuous functions.
%   S = inner_gradP1_gradP1_L2(MESH, PHI, PSI) returns the inner product
%   (grad PHI, grad PSI)_L2.
%
%   V = inner_gradP1_gradP1_L2(MESH, PHI) returns the vector V of grad PHI
%   tested with all gradients of hat functions 
%   V(i) == inner_gradP1_gradP1_L2(MESH, grad PHI, grad HAT(i)).
%
%   M = inner_gradP1_gradP1_L2(MESH) or 
%   M = inner_gradP1_gradP1_L2(MESH, [], [], DIMPHI)
%   returns the matrix M(i,j) = inner_gradP1_gradP1_L2(MESH, HAT(i), HAT(j))
%   [The so called 'stiffness matrix' for P^1 functions].
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~isfield(mesh,'volumes')
    mesh.volumes = getElementVolumes(mesh);
end

if ~exist('phi','var')
    S = inner_gradP1_gradP1_L2_mat(mesh, 1);
elseif ~exist('psi','var')
    S = inner_gradP1_gradP1_L2(mesh, phi, []);
elseif isempty(phi) && isempty(psi)
    S = inner_gradP1_gradP1_L2_mat(mesh, dimPhi);
elseif isempty(phi)
    S = inner_gradP1_gradP1_L2_mat(mesh, 1)*psi;
elseif isempty(psi)
    S = inner_gradP1_gradP1_L2_mat(mesh, 1)*phi;
elseif exist('dimPhi','var') && strcmpi(dimPhi, 'elementwise')
    S = dot(phi, inner_gradP1_gradP1_L2_mat(mesh, 1)*psi, 2);
else
    S = sum(dot(phi, inner_gradP1_gradP1_L2_mat(mesh, 1)*psi, 2));
end
end


function S = inner_gradP1_gradP1_L2_mat(mesh, dim)
% So called 'stiffness matrix'.
if dim==1 && isfield(mesh,'stiffness')
    S = mesh.stiffness;
    return;
end
nC = numCoordinates(mesh);
nE = numElements(mesh);

% (dimMesh+1) hat-functions per Element => nEx(dimMesh+1)x(dimMesh+1) contributions.

hatGrads = getHatGrads(mesh);
[S,I,J] = deal(zeros(nE,dimMesh(mesh)+1,dimMesh(mesh)+1));
for i = 1:dimMesh(mesh)+1
    for j = 1:dimMesh(mesh)+1
		% This works because elementwise P0_P0 scalar product can be interpreted as nE local P0 functions.
		% For higher order elements this is not that simple. See e.g. inner_P1_P1_L2
        S(:,i,j) = inner_P0_P0_L2(mesh,hatGrads{i},hatGrads{j},'elementwise');
        I(:,i,j) = mesh.elements(:,i);
        J(:,i,j) = mesh.elements(:,j);
    end
end
%% Compute the matrix
if dim==1
    S = sparse(I(:),J(:),S(:),nC,nC);
else
    S = sparseBlockDiag(I,J,S,nC,nC,dim);
end
end



%% Alternative (faster, but uglier)
% S = zeros(nE,dimMesh(mesh)+1,dimMesh(mesh)+1);
% for i = 1:dimMesh(mesh)+1
%     for j = 1:dimMesh(mesh)+1
%         if (i<j)
%             [S(:,i,j),S(:,j,i)] = deal(inner_P0_P0_L2(mesh,hatGrads{i},hatGrads{j},'elementwise'));
%         elseif(i==j)
%             S(:,i,i) = inner_P0_P0_L2(mesh,hatGrads{i},hatGrads{j},'elementwise');
%         end
%     end
% end
% idx = ones(dimMesh(mesh)+1,1)*(1:dimMesh(mesh)+1);
% I = elements(:,idx);
% J = elements(:,idx');