function [x, mesh, indicators] = adaptivePoisson(mesh, f, g, uD, nEmax, ...
                                                 rho, discretizer, ...
                                                 quadDeg, postProcessor)
%ADAPTIVEPOISSON    Solves the poisson problem with mixed Dirichlet-Neumann
%boundary conditions using adaptive mesh refinement.
%   Usage:
%   [x, finalMesh, finalIndicators] = ADAPTIVEPOISSON(mesh, f, g, uD, ...
%                                             nEmax, rho, discretizer, ...
%                                             quadDeg, postProcessor)
%
%   postProcessor is an optional function handle with the signature
%   @(mesh, x, indicators) which can be used for plotting or saving the
%   solution. E.g.
%       postProcessor = @surfSolution;
%       postProcessor = @(mesh, x, ind_) saveSolution('example2D', mesh, x);
%
%
%   See also:
%   SOLVEPOISSON
%   SURFSOLUTION
%   SAVESOLUTION
%
%   Author: Josef Kemetmueller - 28.03.2013
dim = @(f) size(f(mesh.coordinates(1,:)),2);
assert(dim(f)==1, 'dim f != 1');
assert(dim(g)==1, 'dim g != 1');
assert(dim(uD)==1, 'dim uD != 1');
while 1
    % Compute discrete solution
    x = solvePoisson(mesh,f,g,uD,discretizer,quadDeg);
    % Compute the refinement indicators
    indicators = estimatorPoissonResidual(mesh, x, f, g, quadDeg);
    % Process solution
    if exist('postProcessor','var')
        postProcessor(mesh, x, indicators);
    end
    % Stopping criterion
    if numElements(mesh) >= nEmax
        break
    end
    % Choose elements for refinement
    [indicators,idx] = sort(indicators,'descend');
    sumeta = cumsum(indicators);
    ell = find(sumeta>=sumeta(end)*rho,1);
    marked = idx(1:ell);
    % Refine Triangulation
    mesh = bisectionRefine(mesh,marked);
end