function saveSolution(name, mesh, x, varargin) %#ok<INUSL>
%SAVESOLUTION    Saves the solution of a simulation
%   Use this function as the postProcessor argument in the adaptivePoisson
%   function.
%   Example:
%      postProcessor = @(mesh, x, ind_) saveSolution('example2D', mesh, x);
%
%   See also:
%   adaptivePoisson
%
% Author: Josef Kemetmueller - 20.03.14

save(sprintf('%s_%d_sol',name, numElements(mesh)), 'x');
save(sprintf('%s_%d_mesh',name, numElements(mesh)), 'mesh');
if nargin>3
    indicators = varargin{1}; %#ok<NASGU>
    save(sprintf('%s_%d_ind',name, numElements(mesh)), 'indicators');
end