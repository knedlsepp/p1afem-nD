function integrals = integrateP0(mesh, P0, varargin)
%INTEGRATEP0    Integrates a piecewise constant function
%   INTEGRAL = INTEGRATEP0(MESH, P0) returns the integral of the function
%   P0 which is given by a nE-by-dimP0 matrix over the given mesh.
%   
%   INTEGRALS = INTEGRATEP0(MESH, P0, 'elementwise') returns a nE-by-dimP0
%   array containing the elementwise integrals.
%
%   INTEGRALS = INTEGRATEP0(MESH, P0, 'patchwise') returns a nC-by-dimP0
%   array containing the patchwise integrals.
%
%   Author: Josef Kemetmueller - 16.12.2013

dimP0 = size(P0,2);
nC = numCoordinates(mesh);
nE = numElements(mesh);
assert(size(P0,1)==nE, ['Integrand must be given as a P^0 function. ',...
                        '[size(elements,1)==size(P0,1)]']);

if isfield(mesh,'volumes')
    volumes = mesh.volumes;
else
    volumes = getElementVolumes(mesh);
end
if (isempty(varargin))
    integrals = volumes'*P0;
else
    switch lower(varargin{1})
        case {'elementwise'}
            integrals = bsxfun(@times, volumes, P0);
        case {'patchwise'}
            integrals = zeros(nC,dimP0);
            for d = 1:dimP0
                for node = 1:dimMesh(mesh)+1
                    integrals(:,d) = integrals(:,d) + ...
                        accumarray(mesh.elements(:,node), ...
                                   bsxfun(@times,volumes,P0(:,d)),[nC,1]);
                end
            end
        otherwise
            error(['Invalid optional argument: ''', ...
                    varargin{1},'''.']);
    end
end
end