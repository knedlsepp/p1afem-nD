function S = inner_P1_P1_L2(mesh, phi, psi, dimPhi)
%INNER_P1_P1_L2    L2-inner product of two elementwise linear, globally
%continuous functions.
%   S = INNER_P1_P1_L2(MESH, PHI, PSI) returns the inner product (PHI,PSI)_L2.
%
%   V = INNER_P1_P1_L2(MESH, PHI) returns the vector V of PHI tested with
%   all hat functions V(i) == INNER_P1_P1_L2(MESH, PHI, HAT(i))
%
%   M = INNER_P1_P1_L2(MESH) or M = INNER_P1_P1_L2(MESH, [], [], DIMPHI)
%   returns the matrix M(i,j) = INNER_P1_P1_L2(MESH, HAT(i), HAT(j))
%   [The so called 'mass matrix' for P^1 functions].
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~isfield(mesh,'volumes')
    mesh.volumes = getElementVolumes(mesh);
end

if ~exist('phi','var')
    S = inner_P1_P1_L2_mat(mesh, 1);
elseif ~exist('psi','var')
    S = inner_P1_P1_L2(mesh, phi, []);
elseif isempty(phi) && isempty(psi)
    S = inner_P1_P1_L2_mat(mesh, dimPhi);
elseif isempty(phi)
    S = inner_P1_P1_L2_mat(mesh, 1)*psi;
elseif isempty(psi)
    S = inner_P1_P1_L2_mat(mesh, 1)*phi;
elseif exist('dimPhi','var') && strcmpi(dimPhi, 'patchwise')
    S = dot(phi, inner_P1_P1_L2_mat(mesh, 1)*psi, 2);
elseif exist('dimPhi','var') && strcmpi(dimPhi, 'elementwise')
    error('TODO');
else
    S = sum(dot(phi, inner_P1_P1_L2_mat(mesh, 1)*psi, 2));
end
end

function S = inner_P1_P1_L2_mat(mesh,dim)
nC = numCoordinates(mesh);
nE = numElements(mesh);

% (dimMesh+1) hat-functions per element => nEx(dimMesh+1)x(dimMesh+1) contributions.
[S,I,J] = deal(zeros(nE,dimMesh(mesh)+1,dimMesh(mesh)+1));

% The following is not the most descriptive formula for the local mass
% matrix, but at least tested for correctness.
inner_HatI_HatJ_ref = factorial(dimMesh(mesh))/factorial(dimMesh(mesh)+2)* ...
                                          (eye(dimMesh(mesh)+1)+ones(dimMesh(mesh)+1));
for i = 1:dimMesh(mesh)+1
    for j = 1:dimMesh(mesh)+1
        S(:,i,j) = mesh.volumes*inner_HatI_HatJ_ref(i,j);
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