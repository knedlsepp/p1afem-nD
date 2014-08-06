function val = norm_P1_L2(mesh, P1)
%NORM_P1_L2    L2-Norm of an elementwise linear, globally continuous function.
%   NORM = NORM_P1_L2(MESH, P1) returns norm of the function P1 induced by
%   the inner_product (.,.)_{L^2} defined by 
%   $(phi, psi)_{L^2} = \int <phi,psi> dx$.
%
%   Author: Josef Kemetmueller - 16.12.2013
val = sqrt(inner_P1_P1_L2(mesh, P1, P1));
end