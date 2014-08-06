function S = inner_P1_P1_h(mesh, phi, psi, dimPhi)
%INNER_P1_P1_H    h-inner product of two piecewise linear globally continous
%functions.
%   S = INNER_P1_P1_H(MESH, PHI, PSI) returns reduced inner_product (.,.)_h
%   defined by $(phi, psi)_h = \int I_h(<phi ,psi>) dx$. ['mass lumping']
%
%   V = INNER_P1_P1_H(MESH, PHI) returns the vector V of PHI tested with
%   all hat functions V(i) == INNER_P1_P1_H(MESH, PHI, HAT(i))
%
%   M = INNER_P1_P1_H(MESH) or M = INNER_P1_P1_H(MESH, [], [], DIMPHI)
%   returns the matrix M(i,j) = INNER_P1_P1_H(MESH, HAT(i), HAT(j))
%   [We could call this the 'reduced mass matrix using mass lumping'].
%
%   Author: Josef Kemetmueller - 16.12.2013

if ~exist('phi','var')
    S = inner_P1_P1_h_mat(mesh, 1);
elseif ~exist('psi','var')
    S = inner_P1_P1_h(mesh, phi, []);
elseif isempty(phi) && isempty(psi)
    S = inner_P1_P1_h_mat(mesh, dimPhi);
elseif isempty(phi)
    % The inner product could be implemented in multiple ways.
    %%% Matrix reuse
    %S = inner_P1_P1_h_mat(mesh, 1)*psi;
    %%% Speed
    S = bsxfun(@times, psi, mesh.betas);
elseif isempty(psi)
    % The inner product could be implemented in multiple ways.
    %%% Matrix reuse
    %S = inner_P1_P1_h_mat(mesh, 1)*phi;
    %%% Speed
    S = bsxfun(@times, phi, mesh.betas);
elseif exist('dimPhi','var') && strcmpi(dimPhi, 'elementwise')
    S = dot(phi, inner_P1_P1_h_mat(mesh, 1)*psi, 2);
else
    % The inner product could be implemented in multiple ways.
    %%% Matrix reuse
    %S = sum(dot(phi, inner_P1_P1_h_mat(mesh, 1)*psi, 2));
    %%% Definition
    S = integrateP1(mesh, dot(phi, psi,2));
end
end

function S = inner_P1_P1_h_mat(mesh,dimPhi)
nC = numCoordinates(mesh);
S = spdiags(repmat(mesh.betas(:),dimPhi,1),0,dimPhi*nC,dimPhi*nC);
end











