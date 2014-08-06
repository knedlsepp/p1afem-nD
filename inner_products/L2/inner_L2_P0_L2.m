function S = inner_L2_P0_L2(mesh, L2, P0, quadDeg, option)
%INNER_L2_P0_L2    Approximate L2-inner product of a L2 function and an
%elementwise constant function.
%   L2 must be given as a function handle, P0 using the elementwise basis.
%
%   V = INNER_L2_P0_L2(MESH, L2) returns the vector V of L2 tested with
%   all constant functions V(i) == INNER_L2_P0_L2(MESH, L2, CHI(i))
%
%   S = INNER_L2_P0_L2(MESH, L2, P0) returns the inner product (L2,P0)_L2.
%
%   S = INNER_L2_P0_L2(MESH, L2, P0, quadDeg) and
%   V = INNER_L2_P0_L2(MESH, L2, [], quadDeg) use quadDeg as quadrature
%   degree.
%
%   Author: Josef Kemetmueller - 16.12.2013

if ~exist('P0','var')
    S = inner_L2_P0_L2_vec(mesh, L2, 3);
elseif isempty(P0)
    S = inner_L2_P0_L2_vec(mesh, L2, quadDeg);
elseif ~exist('quadDeg','var')
    S = inner_L2_P0_L2(mesh, L2, P0, 3);
elseif exist('option','var') && strcmpi(option, 'elementwise')
    S = dot(inner_L2_P0_L2_vec(mesh, L2, quadDeg), P0, 2);
else
    S = sum(dot(inner_L2_P0_L2_vec(mesh, L2, quadDeg), P0, 2));
end
end

function SPs = inner_L2_P0_L2_vec(mesh, L2, quadDeg)
SPs = integrate(mesh, L2, quadDeg, 'elementwise');
end
