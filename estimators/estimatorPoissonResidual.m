function etaR = estimatorPoissonResidual(mesh, x, f, g, quadDeg)
%ESTIMATORPOISSONRESIDUAL    Residual error estimator for the poisson
%equation
%   ETAR = ESTIMATORPOISSONRESIDUAL(mesh, x, f, g, quadDeg) returns the
%   refinement indicators as a P0 function, defined elementwise by
%   eta_T^2 := h_T^2 || f ||_L2(T)^2 
%            + h_T   || J_h(dn U) ||_L2(dT \cap Omega)^2
%            + h_T   || g - dnU ||_L2(dT \cap Gamma_N)^2.
%   h_T being computed via h_T := |T|^(1/dim).
%
%   Works for arbitrary-dimensional meshes.
%
%   Author: Josef Kemetmueller - 28.03.2014

dirichletBd = meshBd(mesh,1);
neumannBd = meshBd(mesh,2);
nE = numElements(mesh);
nD = size(dirichletBd.elements,1);
nB = numBoundaries(mesh);
assert(nB >=2, 'Dirichlet and neumann boundary must be given (possibly empty).');
%%
[face2nodes,element2faces,dirichlet2face,neumann2face] = getHyperfaces(mesh);
% Compute the jumps of the normal derivatives.
Jhdudn = hyperfaceJumpsP0(mesh, gradP1(mesh,x));
skeleton = genMesh(face2nodes, mesh.coordinates);
faceEtaR = inner_P0_P0_L2(skeleton, Jhdudn, Jhdudn, 'elementwise');
% Ignore the contributions on the dirichlet boundary.
faceEtaR(dirichlet2face) = zeros(nD~=0,1); 
% For neumann-faces Jhdudn is the normal derivative, so we don't need to
% recompute the normal derivative.
faceEtaR(neumann2face) = ...
                inner_L2_L2_L2(neumannBd, g, g, quadDeg, 'elementwise') ...
             -2*inner_L2_P0_L2(neumannBd, g, Jhdudn(neumann2face), quadDeg,'elementwise') ...
               +inner_P0_P0_L2(neumannBd, Jhdudn(neumann2face), Jhdudn(neumann2face), 'elementwise');

fintegrals = integrate(mesh, @(X) sum(f(X).^2,2), quadDeg, 'elementwise');
% Compute the residual error estimator.
hT = getElementVolumes(mesh).^(1/dimMesh(mesh));
etaR = hT.^2.*fintegrals + ...
	   hT.*sum(reshape(faceEtaR(element2faces),nE,[]),2);

%% TODO: Dirichlet oscillations: 
%% |T|^(1/d)*(|| \ograd g||_{L^2(E)}^2 - || \Pi_E \ograd g||_{L^2(E)}^2)
% % OLD BLABLA
%    for j = 1:nD
%     face = dirichlet2face(j);
% 
%     %* Oberflaechengradient:
%     coord = mesh.coordinates(face2nodes(face,:),:);
%     [grad,W,sizeF] = ograd(coord,g,quadDeg+1);
%     % Also:  |T|^(1/d)*(|| \ograd g||_{L^2(E)}^2 - || \Pi_E \ograd g||_{L^2(E)}^2)
%     osc1 = W*sum(grad.^2,2);
%     osc2 = sizeF*sum(((W*grad)/sizeF).^2,2);% (W*grad)/sizeF ist das Integralmittel
%     faceEtaR(face) = (1/sizeF)*(osc1 - osc2);
%     % hier durch sizeF dividieren, weiter unten mit sizeT multiplizieren,
%     % dann passts wieder.
% 
%    end 
%% Maybe somewhere along the lines:
% inner_L2_L2_L2(dirichletBd, gradP1....
% Maybe we should implement inner_gradL2_P1_L2 ????