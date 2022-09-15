function [A,p] = getAp(n)
p_src = @(x) 0.01 * (cos(2*pi*x + pi) + 1).* (x<=3);
a = 1;
D = 0.1;

x = linspace(0,10,n+1)';
x_int = x(2:end-1);
p = p_src(x_int);

A = full(gallery('tridiag', n-1, -a/2 - D, 2*D, a/2 - D));