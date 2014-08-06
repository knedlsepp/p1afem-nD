function betas = getBetas(mesh)
%GETBETAS    Factors for the computation of the reduced (mass lumped) inner
%product (.,.)_h.
%   Usage:
%   BETAS = GETBETAS(MESH)
%   
%   It holds BETAS(i)==INTEGRATEP1(mesh, HAT(i)) and using the vector BETAS
%   the reduced inner-product can be computed by
%   $(PHI,PSI)_h = \sum_i betas(i) <PHI(i), PSI(i)>$
%
%   Author: Josef Kemetmueller - 16.12.2013

% The elementwise integral of one hat function is the element volume
% divided by the number of nodes in the element. (For the elements adjacent
% to the node; for the others it is zero)
betas = 1/(dimMesh(mesh)+1)*getPatchVolumes(mesh);
end