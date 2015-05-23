function evaluation = evalPointSources(mesh, pointSources, pointSourceValues)
evaluation = zeros(numCoordinates(mesh),1);
T = triangulation(mesh.elements, mesh.coordinates);

for i = 1:size(pointSources,1)
    barys = T.cartesianToBarycentric((1:numElements(mesh)).', ...
                                     repmat(pointSources(i,:),numElements(mesh),1));
    isInside = all(barys>=0,2);
    evaluation = evaluation + ...
                 accumarray(reshape(T.ConnectivityList(isInside,:),[],1), ...
                            reshape(barys(isInside,:),[],1)*pointSourceValues(i), ...
                            [numCoordinates(mesh),1])/nnz(isInside);
end


