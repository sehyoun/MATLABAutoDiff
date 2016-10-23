function output = dertransp(A,nrow)
%% Fall back file for dertransp
% If you are doing matrix-matrix multiplication, compiling the mex files is
%    highly recommended.
%
% Inputs: dA/dx = derivative of matrix A stacked column-wise
%         nrow = number of rows of matrix A
%
% Output: dA'/dx = derivative of matrix A' stacked column-wise
%
% by SeHyoun Ahn, Oct 2016

[m,l] = size(A);
m = m / nrow;
[i,j,v] = find(A);
aux1 = mod(i,nrow);
aux1(aux1==0) = nrow;
output = sparse(((aux1-1)*m+ceil(i/nrow)),j,v,m*nrow,l);
