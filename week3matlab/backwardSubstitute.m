function [x] = backwardSubstitute(U, y)
% backwardSubstitute.m implements solve of Ux = y for upper triangular U 
% via backward substitution algorithm in MATLAB. Use for demonstration only.
%
% Inputs: 
% U         matrix part of problem
% y         vector part of problem
%
% Outputs:
% x         solution to Ux = y

assert(size(U,1) == size(U,2), 'U is not a square matrix')
assert(size(U,2) == size(y,1), 'U and y dont match')

n = size(U,1);
x = zeros(n,1);

x(n) = y(n) / U(n,n);

for i = (n-1):-1:1
    s = 0;
    for j = (i+1):n
        s = s + U(i,j)*x(j);
    end
    x(i) = (y(i) - s) / U(i,i);
end

end