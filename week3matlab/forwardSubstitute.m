function [y] = forwardSubstitute(L, b)
% forwardSubstitute.m implements solve of Lx = b for lower triangular L via
% forward substitution algorithm in MATLAB. Use for demonstration only.
%
% Inputs: 
% L         matrix part of problem
% b         vector part of problem
%
% Outputs:
% x         solution to Lx = b

assert(size(L,1) == size(L,2), 'L is not a square matrix')
assert(size(L,2) == size(b,1), 'L and b dont match')

n = size(L,1);
y = zeros(n,1);

y(1) = b(1) / L(1,1);

for i = 2:n
    s = 0;
    for j = 1:(i-1)
        s = s + L(i,j)*y(j);
    end
    y(i) = (b(i) - s) / L(i,i);
end

end