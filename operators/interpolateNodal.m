function P1 = interpolateNodal(mesh, f, varargin)
%INTERPOLATENODAL    Nodal interpolant of a function handle.
%   P1 = INTERPOLATENODAL(MESH, f) returns the elementwise linear,
%   globally continuous function P1, whose nodal-values are the evaluations
%   of the function f. Only the nodes used by mesh.elements will be
%   evaluated. So this interpolant can also be used with subsets of the
%   mesh, like the boundaries, e.g.
%   P1 = INTERPOLATENODAL(MESHBD(mesh,1))
%   
%   Works for arbitrary-dimensional meshes.
%
%   See also:
%   INTERPOLATECLEMENT
%   L2PROJECTL2TOP1
%
%   Author: Josef Kemetmueller - 16.12.2013

% We use varargin to discard extra parameters, so it has the same
% signature as interpolateClement, L2ProjectL2ToP1. [hack]

nC = numCoordinates(mesh);
dimf = size(f(mesh.coordinates(1,:)),2);

nodes = unique(mesh.elements); 
P1 = zeros(nC,dimf);
P1(nodes,:) = f(mesh.coordinates(nodes,:));