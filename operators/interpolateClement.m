function P1 = interpolateClement(mesh, f, quadDeg)
%INTERPOLATECLEMENT    Clement interpolation of a function handle.
%   P1 = INTERPOLATECLEMENT(MESH, f) returns the elementwise linear,
%   globally continuous function P1, whose nodal-values are the mean values
%   of the corresponding node-patches.
%
%   See also:
%   L2PROJECTL2TOP1
%   INTERPOLATENODAL
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~exist('quadDeg','var');
    quadDeg = 3;
end
P1 = bsxfun(@rdivide, integrate(mesh,f,quadDeg,'patchwise'), getPatchVolumes(mesh));
end
