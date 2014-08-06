function S = inner_L2_P1_L2(mesh, L2, P1, quadDeg, option)
%INNER_L2_P1_L2    Approximate L2-inner product of a L2 function and an
%elementwise linear, globally continuous function.
%   L2 must be given as a function handle, P1 using the nodal basis.
%
%   V = INNER_L2_P1_L2(MESH, L2) returns the vector V of L2 tested with
%   all hat functions V(i) == INNER_L2_P1_L2(MESH, L2, HAT(i))
%
%   S = INNER_L2_P1_L2(MESH, L2, P1) returns the inner product (L2,P1)_L2.
%
%   S = INNER_L2_P1_L2(MESH, L2, P1, quadDeg) and
%   V = INNER_L2_P1_L2(MESH, L2, [], quadDeg) use quadDeg as quadrature
%   degree.
%
%   Author: Josef Kemetmueller - 16.12.2013
if ~isfield(mesh,'volumes')
    mesh.volumes = getElementVolumes(mesh);
end

if exist('option','var') && strcmpi(option, 'elementwise')
    S = inner_L2_P1_L2_vec_elementwise(mesh, L2, P1, quadDeg);
elseif exist('option','var') && strcmpi(option, 'patchwise')
    S = dot(inner_L2_P1_L2_vec(mesh, L2, quadDeg), P1, 2);
elseif ~exist('P1','var')
    S = inner_L2_P1_L2_vec(mesh, L2, 3);
elseif isempty(P1)
    S = inner_L2_P1_L2_vec(mesh, L2, quadDeg);
elseif ~exist('quadDeg','var')
    S = inner_L2_P1_L2(mesh, L2, P1, 3);
else
    S = sum(dot(inner_L2_P1_L2_vec(mesh, L2, quadDeg), P1, 2));
end
end

function SPs = inner_L2_P1_L2_vec(mesh, L2, quadDeg)

dimL2 = size(L2(mesh.coordinates(1,:)),2);
nC = numCoordinates(mesh);


[W, V] = simplexQuadratureRule(dimMesh(mesh), quadDeg);
nQ = size(V,1);
SPs = zeros(nC,dimL2);

for j = 1:nQ
    XYZ = barycentricToCartesians(mesh,V(j,:));
    L2XYZ = L2(XYZ);
    for d = 1:dimL2
        for node = 1:dimMesh(mesh)+1
            SPs(:,d) = SPs(:,d) + accumarray(mesh.elements(:,node), ...
                                W(j)*V(j,node)*mesh.volumes.*L2XYZ(:,d),[nC,1]);
        end
    end
end

end

function SPs = inner_L2_P1_L2_vec_elementwise(mesh, L2, P1, quadDeg)

dimL2 = size(L2(mesh.coordinates(1,:)),2);
nE = numElements(mesh);

[W, V] = simplexQuadratureRule(dimMesh(mesh), quadDeg);
nQ = size(V,1);
SPs = zeros(nE,1);

for j = 1:nQ
    XYZ = barycentricToCartesians(mesh,V(j,:));
    L2XYZ = L2(XYZ);
    for d = 1:dimL2
        for node = 1:dimMesh(mesh)+1
            SPs(:) = SPs(:) + W(j)*V(j,node)...
                              *mesh.volumes.*L2XYZ(:,d).*P1(mesh.elements(:,node),d);
        end
    end
end

end
