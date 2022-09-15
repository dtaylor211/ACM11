function [xn] = fpi(g,x0,tau,Nmax)
% Solves x = g(x) via fixed point iteration
% Inputs: 
%    g     function to find a fixed point of
%    x0    initial guess
%    tau   tolerance for convergence
%    Nmax  max number of iterations
% Outputs:
%    xn    sequence of iterates

xn = zeros(length(x0),Nmax);    % pre allocate space for iterates
xn(:,1) = g(x0);                % initialize

for i = 2:Nmax
    xn(:,i) = g(xn(:,i-1));
    if norm(xn(:,i)-xn(:,i-1))/norm(xn(:,i-1)) < tau
        break
    end
end
xn = xn(:,1:i); % chop off extra zeros
end