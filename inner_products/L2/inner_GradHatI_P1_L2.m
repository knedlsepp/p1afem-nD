function P1TGradHat = inner_GradHatI_P1_L2(mesh, P1)
%INNER_GRADHATI_P1_L2    Elementwise linear function tested by hat gradients.
%   V = INNER_GRADHATI_P1_L2(mesh, P1) returns the vector of P1 tested by
%   all gradients of hat functions
%   V(i) == (P1, grad hat(i))_L2
%
%   Author: Josef Kemetmueller - 16.12.2013


%TODO: Describe P1 multidimensional.
%TODO: inner_gradP1_P1_L2-matrix instead of this. But keep in mind, that LLG uses this code. so if you change this file. then change the llg-code too!
nC = numCoordinates(mesh);
P1TGradHat = zeros(nC,1);
hatGrads = getHatGrads(mesh);
for node = 1:dimMesh(mesh)+1
    P1TGradHat = P1TGradHat + accumarray(mesh.elements(:,node), ...
        inner_P0_P1_L2(mesh, hatGrads{node}, P1, 'elementwise'), [nC,1]);
end
