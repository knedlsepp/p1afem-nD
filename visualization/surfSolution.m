function surfSolution(mesh, solution, varargin)
%SURFSOLUTION    Plots the solution of a simulation
%   Use this function as the postProcessor argument in the adaptivePoisson
%   function.
%   Example:
%      postProcessor = @surfSolution;
%
%   See also:
%   adaptivePoisson
%
% Author: Josef Kemetmueller - 20.03.14

figure(1); clf;
%subplot(1,nargin,1);
surfMesh(mesh, solution);
title(sprintf('Solution. #elements: %d',numElements(mesh)));
colorbar;
%subplot(1,nargin,2);
figure(2);
surfMesh(mesh);
title(sprintf('Boundary mesh. #elements: %d',numElements(mesh)));
if nargin>2
    for j = 1:nargin-2
%        subplot(1,nargin,j+2);
        figure(3);
        surfMesh(mesh, interpolateClementP0(mesh,varargin{j}));
        title(sprintf('Error estimator. #elements: %d',numElements(mesh)));
        colorbar;
        drawnow;
    end
end

end