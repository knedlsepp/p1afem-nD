function val = norm_L2_L2(mesh, L2)
%NORM_L2_L2    L2-Norm of a function handle
%   NORM = NORM_L2_L2(MESH, f) returns norm of the function f induced by
%   the inner_product (.,.)_{L^2} defined by 
%   $(f, f)_{L^2} = \int <f,f> dx$.
%
%   Author: Josef Kemetmueller - 16.12.2013
val = sqrt(inner_L2_L2_L2(mesh, L2, L2));
end