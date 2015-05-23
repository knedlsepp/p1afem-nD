function P1 = L2ProjectL2ToP1(mesh, L2, quadDeg)
%L2PROJECTL2TOP1 L2-Projection from L2 to P1.
%   P1 = L2PROJECTL2TOP1(MESH, L2) returns the L2-projection onto the space
%   of elementwise linear, globally continuous functions.
%
%   See also:
%   L2PROJECTL2TOP0
%   INTERPOLATECLEMENT
%   INTERPOLATENODAL
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~exist('quadDeg','var');
    quadDeg = 3;
end
%% We could do this directly if we knew the mesh was clean.
%P1 = inner_P1_P1_L2(mesh)\inner_L2_P1_L2(mesh, L2, [], quadDeg);

%% To be able to project onto boundaries or unclean meshes:
% If not all coordinates are used by the mesh we clean the mesh
% before, otherwise the matrix will be singular. (It still seems to work
% without cleaning though. Due to the superb sparse solver...)
nC = numCoordinates(mesh);
dimL2 = size(L2(mesh.coordinates(1,:)),2);
[cleanedmesh, ~, new2old] = cleanMesh(mesh);
P1 = zeros(nC,dimL2);
if ~isempty(cleanedmesh.elements)
    P1(new2old,:) = inner_P1_P1_L2(cleanedmesh)\inner_L2_P1_L2(cleanedmesh, L2, [], quadDeg);
end