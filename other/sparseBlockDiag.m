function S = sparseBlockDiag(I,J,S,m,n,repnum)
%SPARSEBLOCKDIAG Sparse block diagonal matrix.
%   S = SPARSEBLOCKDIAG(I,J,S,m,n,repnum) produces a block diagonal matrix,
%   with repnum diagonal blocks consisting of the matrix SPARSE(I,J,S,m,n)
%
%   Author: Josef Kemetmueller - 16.12.2013
ILocal2Global = @(d,nodes) bsxfun(@plus,m*(d-1),reshape(nodes,[],1));
JLocal2Global = @(d,nodes) bsxfun(@plus,n*(d-1),reshape(nodes,[],1));

S = repmat(S(:),1,repnum);
I = ILocal2Global(1:repnum,I);
J = JLocal2Global(1:repnum,J);

S = sparse(I,J,S,m*repnum,n*repnum);

end