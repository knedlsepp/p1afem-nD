
function val = norm_P1_h(mesh, P1)
%NORM_P1_H    h-Norm of an elementwise linear, globally continuous function.
%   NORM = NORM_P1_H(MESH, P1) returns the norm of the function P1 induced
%   by the reduced inner_product (.,.)_h defined by
%   $(phi, psi)_h = \int I_h(<phi,psi>) dx$.
%
%   Author: Josef Kemetmueller - 16.12.2013
val = sqrt(inner_P1_P1_h(mesh, P1, P1));
end