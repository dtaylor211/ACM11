function [x] = solveLU(A, b)
% solveLU.m implements linear solve of Ax = b via LU decomposition and
% subsequent forward/backward substitution. Use of demonstration only.
%
% Inputs: 
% A         matrix part of problem
% b         vector part of problem
%
% Outputs:
% x         solution to Ax = b

assert(size(A,1) == size(A,2), 'A is not a square matrix')
assert(size(A,2) == size(b,1), 'A and b dont match')

[L, U] = getLU(A);
y = forwardSubstitute(L, b);
x = backwardSubstitute(U,y);

end