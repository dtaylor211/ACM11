function [x] = solveCholesky(A, b)
% solveCholesky.m implements linear solve of Ax = b via Cholesky 
% factorization in MATLAB. Use for demonstration only.
%
% Inputs: 
% A         matrix part of problem
% b         vector part of problem
%
% Outputs:
% x         solution to Ax = b

assert(size(A,1) == size(A,2), 'A is not a square matrix')
assert(size(A,2) == size(b,1), 'A and b dont match')

L = cholesky(A);
y = forwardSubstitute(L, b);
x = backwardSubstitute(L',y);

end