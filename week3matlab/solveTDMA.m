function [x] = solveTDMA(A, d)
% solveTDMA.m implements linear solve of Ax = d for tridiagonal A via
% Thomas algorithm (aka TriDiagonal Matrix Algorithm).
% 
% Inputs: 
% A         matrix part of problem
% d         vector part of problem
%
% Outputs:
% x         solution to Ax = d

assert(size(A,1) == size(A,2), 'A is not a square matrix')
assert(size(A,2) == size(d,1), 'A and d dont match')
 
n = size(A,1);
x = zeros(n,1);
e = zeros(n,1);
f = zeros(n,1);

b = A(1,1);
c = A(1,2);

e(1) = c / b;
f(1) = d(1) / b;

for i = 2:(n-1)
    a = A(i,i-1);
    b = A(i,i);
    c = A(i,i+1);
    e(i) = c / (b - a*e(i-1));
end

for i = 2:n
    a = A(i,i-1);
    b = A(i,i);
    f(i) = (d(i) - a*f(i-1)) / (b - a*e(i-1));
end

x(n) = f(n);
for i = (n-1):-1:1
    x(i) = f(i) - e(i)*x(i+1);
end

end

