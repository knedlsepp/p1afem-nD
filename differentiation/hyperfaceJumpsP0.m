function jumps = hyperfaceJumpsP0(mesh, P0)
%HYPERFACEJUMPSP0    Face jumps of the normal component of an elementwise
%constant function
%   JUMPS = HYPERFACEJUMPSP0(MESH, P0) returns an array JUMPS corresponding
%   to the faces returned by GETHYPERFACES. P0 is a dimSpace-vector-valued
%   elementwise constant function. JUMPS is a P0 function defined on its
%   hyperfaces("skeleton"). For boundary faces it yields the outwards jump.
%   JUMPS is NOT scaled by the areas of the hyperfaces. As
%   the skeleton is also a simplicial complex you can use the
%   integration/inner_product routines to integrate the jumps, to get the
%   the jumps scaled with their corresponding areas.
%     Example: Jumps of normal derivative of a P1 function in L1/L2 norm.
%
%       Jhdudn = hyperfaceJumpsP0(mesh, gradP1(mesh,x));
%       skeleton = genMesh(face2nodes, mesh.coordinates);
%       N1 = integrateP0(skeleton, Jhdudn, 'elementwise');
%       N2 = integrateP0(skeleton, Jhdudn.^2, 'elementwise');
%       N2 = inner_P0_P0_L2(skeleton, Jhdudn, Jhdudn, 'elementwise');
%
%   Works for arbitrary-dimensional meshes.
%
%   See also:
%	GETHYPERFACES
%
%   Author: Josef Kemetmueller - 16.12.2013

[face2nodes,element2faces] = getHyperfaces(mesh); % Faces not oriented
nF = size(face2nodes,1);
% Correct orientation of normals at boundary follows from getNormals
normals = getNormals(mesh, 'outwards');
% Normal jumps
jump = cellfun(@(N) dot(N,P0,2), normals, 'UniformOutput', false);
jumps = accumarray(element2faces(:),cell2mat(jump'),[nF,1]);