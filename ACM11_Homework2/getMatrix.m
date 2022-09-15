function A = getMatrix(N, L)
    dx = L/N;
    d1 = diag(1/dx^2 + zeros(N-2,1),-1);
    d2 = diag(-2/dx^2 + zeros(N-1,1));
    d3 = diag(1/dx^2+zeros(N-2,1),1);
    A = d1+d2+d3;
end