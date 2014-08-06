function [W, V] = simplexQuadratureRule(dim, n)
%SIMPLEXQUADRATURERULE    Quadrature points for simplices.
%   [WEIGHTS, POINTS] = SIMPLEXQUADRATURERULE(DIM, QUADDEG) returns
%   quadrature weights and points in barycentric coordinates for a simplex
%   of dimension DIM using QUADDEG^DIM quadrature points.
%   WEIGHTS is a row vector of quadrature weights corresponding to the rows
%   of the matrix POINTS. Which are the quadrature points given in
%   barycentric coordinates. 
%   If a function f can evaluate multiple points at once you can simply use
%       integral = volume*W*f(V)
%   meaning $\int f dx = volume \sum_i w_i f(V_i)$.
%
%   Otherwise use the following code 
%       integral = 0;
%       for i = 1:length(W)
%           integral = integral + volume*W(i)*f(V(i,:))
%       end
%   The quadrature rule is exact for polynomials of degree 2*QUADDEG-DIM.
%
%   See also:
%	BARYCENTRICTOCARTESIANS
%
%   Author: Josef Kemetmueller - 16.12.2013

% One-point quadrature. 
if (isempty(n) || (n == 1))
    V = 1/(dim+1)*ones(1,dim+1);
    W = 1;
    return;
end
% We use tensor quadrature based on 1D-gaussian quadrature and
% integration by substitution via the affine transformation TRAFO and its
% determinant TRAFODET.
[X_1D,W_1D] = gauss1D(n,0,1);
switch dim
    case -1 %empty mesh
        W = []; V = []; return;
    case 0
        W = 1; V = 1;
    case 1
        W = W_1D;
        V = [X_1D, 1-X_1D];
    case 2
        trafo = @(x,y) [x,(1-x).*y];
        trafodet = @(x,y) (1-x);
        % Tensor quadrature
        X = repmat(X_1D',[n,1]);
        W = repmat(W_1D,[n,1]);
        Y = X';
        W = W.*W';
        % Affine transformation of quadrature rule
        XY = trafo(X(:),Y(:));
        V = [XY,1-sum(XY,2)];
        W = W.*abs(trafodet(X,Y));        
    case 3
        trafo = @(x,y,z) [x,(1-x).*y,(1-x).*(1-y).*z];
        trafodet = @(x,y,z) (1-x).*(1-x).*(1-y);
        % Alternative for edge-singularities
        % trafo = @(x,y,z) [x.*(1-y), x.*y, (1-x).*z];
        % trafodet = @(x,y,z) (1-y).*x.*(1-x);
        % Tensor quadrature
        X = repmat(X_1D',[n,1,n]);
        W = repmat(W_1D,[n,1,n]);
        Y = shiftdim(X,1);
        Z = shiftdim(X,2);
        W = W.*shiftdim(W,1).*shiftdim(W,2);
        % Affine transformation of quadrature rule
        XYZ = trafo(X(:),Y(:),Z(:));
        V = [XYZ,1-sum(XYZ,2)];
        W = W.*abs(trafodet(X,Y,Z));
    otherwise
        W0 = repmat(W_1D,[n,1,n+zeros(1,dim-2)]);
        W = W0;
        for j = 1:dim-1
            W = W.*shiftdim(W0,j);
        end
        X = cell(1,dim);
        X0 = repmat(X_1D',[n,1,n+zeros(1,dim-2)]);
        for d1 = 1:dim
            X{d1} = shiftdim(X0,d1-1);
            for d2 = 1:d1-1
                X{d1} = X{d1}.*(1-shiftdim(X0,d2-1));
                W = W.*abs(1-shiftdim(X0,d2-1));
            end
        end
        XYZ = cell2mat(cellfun(@(X) reshape(X,[],1), X, 'UniformOutput',false));
        V = [XYZ,1-sum(XYZ,2)];
end
% As the above yields a quadrature rule for the n-dimensional simplex, it
% holds sum(W) = 1/factorial(dim) [= volumeOf(nDimSimplex)]. But as we want
% to use the function via integral = volume*W*f(V) we have to scale so
% that sum(W) == 1.
W = reshape(factorial(dim)*W,1,[]);
end

function [nodes,weights] = gauss1D(n,varargin)
beta = (1:n-1)./sqrt((2*(1:n-1)).^2-1);
A = diag(beta,-1)+diag(beta,1);
[eigenvector,nodes] = eig(A);
[nodes,idx] = sort(diag(nodes));
weights = 2*eigenvector(1,idx).^2; 

if nargin >= 3
    a = varargin{1};
    b = varargin{2};
    weights = 0.5*abs(b-a)*weights;
    nodes = 0.5*( a+b + nodes*(b-a) );
end
end
