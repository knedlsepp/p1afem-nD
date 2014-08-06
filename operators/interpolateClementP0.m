function P1 = interpolateClementP0(mesh,P0)
%INTERPOLATECLEMENTP0 Clement interpolation of a P0 function.
%   P1 = INTERPOLATECLEMENTP0(MESH, P0) returns the elementwise linear,
%   globally continuous function P1, whose nodal-values are the mean values
%   of the corresponding node-patches.
%
%   Author: Josef Kemetmueller - 16.12.2013

% Integral mean of the patch
P1 = bsxfun(@rdivide, integrateP0(mesh,P0,'patchwise'), getPatchVolumes(mesh));
end