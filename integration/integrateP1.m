function integrals = integrateP1(mesh, P1, varargin)
%INTEGRATEP1    Integrates a piecewise linear function
%   INTEGRAL = INTEGRATEP1(MESH, P1) returns the integral of the function
%   P1 which is given by a nC-by-dimP1 matrix over the given mesh,
%   which is a structure containing the arrays ELEMENTS and COORDINATES.
%   These represent the mesh by standard simplex-vertex format.
%   
%   INTEGRALS = INTEGRATEP1(MESH, P1, 'elementwise') returns a nE-by-dimP1
%   array containing the elementwise integrals.
%
%   INTEGRALS = INTEGRATEP1(mesh, P1, 'patchwise') returns a nC-by-dimP1
%   array containing the patchwise integrals.
%
%   Author: Josef Kemetmueller - 16.12.2013
nC = numCoordinates(mesh);
assert(size(P1,1)==nC, ['Integrand must be given as a P^1 function. ', ...
                        '[size(coordinates,1) == size(P1,1)]']);
%% Reuse the available functions
integrals = integrateP0(mesh, L2ProjectP1ToP0(mesh, P1), varargin{:});