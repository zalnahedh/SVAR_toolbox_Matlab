function output = logAbsDet(X)
%------------------------------------------------------------------------
% INPUT:  square matrix X
% OUTPUT: real number that is equal to logarithm of absolute value of
%         determinant
% Function uses LU decomposition:
%   Any permutation matrix has determinant ±1, depending on the parity of the permutation
%   To find the determinant of an upper triangular or lower triangular matrix, take the product of the diagonal entries.
%   If A=PLU, then det(A)=det(P)det(L)det(U)
%   Since L has ones as diagonal entries and we take absolute value, only U
%   is needed
%-------------------------------------------------------------------------

[~,U,~]=lu(X);
n=size(U,1);
x=0.0;
for i=1:n
    if U(i,i) == 0.0
        x=-inf;
        return;
    end
    x=x+log(abs(U(i,i)));
end
output=x;
end
    