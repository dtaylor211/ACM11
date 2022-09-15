% 16.S685 getLU.m PROVIDED CODE
% calculates the LU-decomposition of A without pivoting
%
% Inputs: 
% A         matrix to perform LU-decomposition on
%
% Outputs:
% L         lower triangular matrix
% U         upper triangular matrix

function [L,U] = getLU(A)
% Ensures A is square
assert(size(A,1) == size(A,2), 'A is not a square matrix')
 
n = size(A,1);
L = eye(n);
U = zeros(n,n);

for j = 1:n
    for i = 1:j
        s1 = 0;
        for k = 1:i
            s1 = s1 + U(k,j)*L(i,k);
        end
        U(i,j) = A(i,j) - s1;
    end
    for i = (j+1):n
        s2 = 0;
        for k = 1:(j-1)
            s2 = s2 + U(k,j)*L(i,k);
        end
        L(i,j) = (A(i,j)-s2) / U(j,j);
    end
end

end

