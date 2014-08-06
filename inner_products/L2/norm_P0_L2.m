function val = norm_P0_L2(mesh, P0)
%NORM_P0_L2    L2-Norm of an elementwise constant function.
%   NORM = NORM_P0_L2(MESH, P1) returns norm of the function P1 induced by
%   the inner_product (.,.)_{L^2} defined by 
%   $(CHI_1, CHI_2)_{L^2} = \int <CHI_1,CHI_2> dx$.
%
%   Author: Josef Kemetmueller - 16.12.2013
val = sqrt(inner_P0_P0_L2(mesh, P0, P0));
end