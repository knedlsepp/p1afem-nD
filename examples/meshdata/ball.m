function mesh = ball(numIt)
%BALL    Generates the approximate mesh of a ball.
%   BALL(NUMIT) returns an approximate mesh of the unit ball bisecting
%   the initial mesh numIt times.
%
%   Author: Josef Kemetmueller - 16.12.2013

%%
bisFactor = 100/100; % How much to bisect each time.
startTriangulation = 1;
switch startTriangulation
    case 1
    [X,Y,Z] = meshgrid(-1:2:1,-1:2:1,-1:2:1);
    coordinates = [X(:), Y(:), Z(:)];
    case 2
    coordinates = ...
    [ 0, 0, 0;
     -1, 0, 0;
      1, 0, 0;
      0, 1, 0;
      0,-1, 0;
      0, 0, 1;
      0, 0,-1];
    case 3
    coordinates = ...
    [0.000, 0.000, 1.000 
     0.943, 0.000,-0.333 
    -0.471, 0.816,-0.333 
    -0.471,-0.816,-0.333 ];
end
mesh = genBisectionMesh(genMesh(delaunay(coordinates), coordinates,'any'));
lengths = @(X) sqrt(sum(X.^2,2));
normalize = @(X) bsxfun(@rdivide, X, lengths(X));


for i = 1:numIt
    nE = numElements(mesh);
    volumes = getElementVolumes(mesh);
    [~, I_sort] = sort(volumes);
    marked = I_sort(fix(nE*bisFactor):end);
    mesh = bisectionRefine(mesh, marked);

    boundary = getBoundary(mesh);
    bNodes = unique(boundary(:));
    mesh.coordinates(bNodes,:) = normalize(mesh.coordinates(bNodes,:));
end
end