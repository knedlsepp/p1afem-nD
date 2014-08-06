function [integrals, varargout] = integrate(mesh, f, quadDeg, varargin)
%INTEGRATE    Integrates a general function given as a function handle.
%   INTEGRAL = INTEGRATE(MESH, f, quadDeg) returns the integral of the
%   function f which is given by a function handle over the given mesh.
%   quadDeg is an integer representing the number of quadrature points used
%   in the underlying tensor gaussian quadrature.
%
%   INTEGRALS = INTEGRATE(MESH, f, quadDeg, 'elementwise') returns a
%   nE-by-dimf array containing the elementwise integrals.
%
%   INTEGRALS = INTEGRATE(mesh, f, quadDeg, 'patchwise') returns
%   nC-by-dimf array containing the patchwise integrals.
%
%   WARNING: The function f must be given in a manner, that n points can
%   be evaluated simultaneously given as a n-by-dimSpace matrix.
%   If you don't have a vectorized version of the function f it can be
%   integrated using the [slower] option 'for':
%   INTEGRALS = INTEGRATE(..., 'for')
%
%   Example 1D:
%       n = 9;
%       elements = [(1:n-1)', (2:n)'];
%       coordinates = linspace(0,pi/2, n)';
%       mesh = genMesh(elements, coordinates);
%       int = integrate(mesh, @sin)
%
%   Example 2D:
%       [X,Y] = meshgrid(0:1,0:1);
%       coordinates = [X(:), Y(:)];
%       elements = delaunay(coordinates);
%       mesh = genMesh(elements, coordinates);
%       int = integrate(mesh, @(X) X(:,1).^2)
%
%   The used quadrature weights W and points C can be returned using
%   [INTEGRAL, W, C] = INTEGRATE(mesh, f, quadDeg)
%   and evaluated using: 
%   integrals = bsxfun(@times, getElementVolumes(mesh), ...
%                      reshape(f(reshape(C,[],size(C,3))),[],length(W)))*W
%
%   Works for arbitrary-dimensional meshes.
%
%   See also:
%   GENMESH
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~exist('quadDeg','var') || isempty(quadDeg)
    quadDeg = 3;
end
dimf = size(f(mesh.coordinates(1,:)),2);
assert(size(f(mesh.coordinates(1,:)),1)==1,'f values must be row vectors.');

if any(strcmpi('for',varargin))
    f = @(X) manualVectorization(f,X,dimf);
else
    assert(size(f([mesh.coordinates(1,:);mesh.coordinates(1,:)]),1) == 2, ...
                            'Function not vectorized. Use ''for'' option.');
end

nC = numCoordinates(mesh);
nE = numElements(mesh);
%%
integral_elements = zeros(nE,dimf);
[W, V] = simplexQuadratureRule(dimMesh(mesh), quadDeg);
if nargout>1
    varargout{1} = W';
end
if nargout>2
    varargout{2} = zeros(nE,length(W),3);
end
volumes = getElementVolumes(mesh);
nQ = length(W);
for j = 1:nQ
    XYZ = barycentricToCartesians(mesh,V(j,:));
    if nargout>2
        varargout{2}(:,j,:) = XYZ;
    end
    fXYZ = f(XYZ);
    for d = 1:dimf
        integral_elements(:,d) = integral_elements(:,d) + W(j)*volumes.*fXYZ(:,d);
    end
end

if (isempty(varargin) || strcmpi(varargin{1},'for'))
    integrals = sum(integral_elements,1);
else
    switch lower(varargin{1})
        case {'elementwise'}
            integrals = integral_elements;
        case {'patchwise'}
            integrals = zeros(nC,dimf);
            for d = 1:dimf
                for node = 1:dimMesh(mesh)+1
                    integrals(:,d) = integrals(:,d) + ...
                                     accumarray(mesh.elements(:,node), ...
                                                integral_elements(:,d),[nC,1]);
                end
            end
        otherwise
            error(['Invalid optional argument: ''', ...
                    varargin{1},'''.']);
    end
end
end

function Y = manualVectorization(f,X,dimf)
% If f is not vectorized, we use a loop to calculate the values.
Y = zeros(size(X,1), dimf);
for i = 1:size(X,1)
    Y(i,:) = f(X(i,:));
end
end
