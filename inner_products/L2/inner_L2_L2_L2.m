function S = inner_L2_L2_L2(mesh, f, g, quadDeg, option)
%INNER_L2_L2_L2    Approximate L2-inner product of two L2 functions.
%   S = INNER_L2_L2_L2(MESH, F, G) returns the L2-scalar product of F and G.
%   The functions must be given as function handles.
% 
%   S = INNER_L2_L2_L2(MESH, F, G, quadDeg) uses quadDeg as quadrature
%   degree.
% 
%   Ss = INNER_L2_L2_L2(MESH, F, G, quadDeg, 'elementwise') yields the
%   elementwise L2-scalar products.
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~exist('quadDeg','var') || isempty(quadDeg)
    quadDeg = 3;
end
  
if exist('option','var')
    S = integrate(mesh, @(X) dot(f(X),g(X),2), quadDeg, option);
else
    S = integrate(mesh, @(X) dot(f(X),g(X),2), quadDeg);
end
end