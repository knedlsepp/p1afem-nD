function surfMesh(mesh,varargin)
%SURFMESH    Plots a mesh data structure
%   surfMesh(mesh) will plot a one, two or three dimensional mesh.
%   surfMesh(mesh, P1) will plot a piecewise linear function P1
%   corresponding to the nodes of the mesh.
%   
%   Author: Josef Kemetmueller - 28.03.2013

assert(dimMesh(mesh)<=3, 'No visualization implemented for dimMesh=%d', ...
                         dimMesh(mesh));
assert(dimSpace(mesh)<=3, 'No visualization implemented for dimSpace=%d', ...
                          dimMesh(mesh));           
was_hold = ishold();
switch dimMesh(mesh)
    case 1
        surfMesh1D(mesh,varargin{:});
    case 2
        surfMesh2D(mesh,varargin{:});
    case 3
        surfMesh3D(mesh,varargin{:});
end
if ~was_hold
    hold off;
end
end
%% 1D
function surfMesh1D(mesh, varargin)

if dimSpace(mesh)==1
        [C,I] = sort(mesh.coordinates);    
        if nargin>1 && isnumeric(varargin{1})
            plot(C, varargin{1}(I));
        else
            plot(C, zeros(size(C)));
        end
else
    nC = numCoordinates(mesh);
    if nargin>1 && isnumeric(varargin{1})
        C = [mesh.coordinates, reshape(varargin{1},[],1)];
        trisurf(mesh.elements(:,[1,2,2]),C(:,1), C(:,2), C(:,3), ...
                'facecol','no','edgecol','interp','linew',3);
    else
        C = [mesh.coordinates, zeros(nC,3-dimSpace(mesh))];
        trisurf(mesh.elements(:,[1,2,2]),C(:,1), C(:,2), C(:,3), ...
               'facecol','no','edgecol','flat','linew',3);
    end

end

end

%% 2D
function surfMesh2D(mesh,varargin)
nC = numCoordinates(mesh);
C = @(i) mesh.coordinates(:,i);
if dimSpace(mesh)==3
    if nargin>1 && isnumeric(varargin{1})
        trisurf(mesh.elements,C(1), C(2), C(3), varargin{1},'facecolor','interp','edgecolor', 'none');
    else
        trisurf(mesh.elements,C(1), C(2), C(3), zeros(nC,1),'facecolor','flat');
    end
else
    if nargin>1 && isnumeric(varargin{1})
        trisurf(mesh.elements,C(1), C(2), varargin{1},'facecolor','interp','edgecolor', 'none');
    else
        trisurf(mesh.elements,C(1), C(2), zeros(nC,1),'facecolor','interp');
    end
end
axis equal;
axis vis3d;
view(15,22);
end
%% 3D
function surfMesh3D(mesh,varargin)
nB = numBoundaries(mesh);
if nB ==0
    nB = 1;
    mesh.bd(1).elements = getBoundary(mesh);
end
colors = colormap;
cnum = linspace(1,size(colors,1), nB);
C = @(i) mesh.coordinates(:,i);
midpoints = @(EL, C) 1/3*(C(EL(:,1),:)+C(EL(:,2),:)+C(EL(:,3),:));
normalvectors = @(EL, C) cross(C(EL(:,2),:)-C(EL(:,1),:),C(EL(:,3),:)-C(EL(:,1),:),2);
if nargin>1 && isnumeric(varargin{1})
    x = varargin{1};
    for j = 1:nB
        trisurf(mesh.bd(j).elements, C(1), C(2), C(3), reshape(x,[],1), 'facecolor', 'interp', 'FaceLighting', 'phong','edgecolor', 'none');
        hold on;
    end
else
    for j = 1:nB
        trisurf(mesh.bd(j).elements, C(1), C(2), C(3), 'FaceColor', colors(cnum(j),:),'faceAlpha', 0.5);
        hold on;
    end
    %tetramesh(mesh.elements, mesh.coordinates, 'FaceColor', 'flat','FaceAlpha', 0.1);
    %tetramesh(mesh.elements, mesh.coordinates, 'FaceColor', 'none', 'EdgeLighting', 'phong', 'EdgeAlpha', 0.3);
end

if nargin>1 && any(cellfun(@(C) strcmpi(C,'quiver'),varargin))
    for j = 1:nB
        MP = midpoints(mesh.bd(j).elements, mesh.coordinates);
        NV = normalvectors(mesh.bd(j).elements, mesh.coordinates);
        quiver3(MP(:,1),MP(:,2),MP(:,3),NV(:,1),NV(:,2),NV(:,3), 'color', 'k');
    end
end
axis vis3d;
view([-214,40]);

end