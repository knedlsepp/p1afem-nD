function P0 = L2ProjectL2ToP0(mesh, L2, quadDeg)
%L2PROJECTL2TOP0 L2-Projection from L2 to P0.
%   P1 = L2PROJECTL2TOP1(MESH, L2) returns the L2-projection onto the space
%   of elementwise constant functions.
%
%   See also:
%   L2PROJECTL2TOP1
%   INTERPOLATECLEMENT
%   INTERPOLATENODAL
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~exist('quadDeg','var');
    quadDeg = 3;
end
if ~isempty(mesh.elements)
    P0 = inner_P0_P0_L2(mesh)\inner_L2_P0_L2(mesh, L2, [], quadDeg);
end