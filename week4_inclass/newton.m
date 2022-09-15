function xn = newton(f,J,x0,tau,Nmax)
% Solves the nonlinear equation given by f(x) = 0 via Newton iteration
% Inputs: 
%    f         function to find root of
%    J         Jacobian of f
%    x0        initial guess
%    tau       convergence tolerance
%    Nmax      maximum number of iterations
% Outputs:
%    xn        matrix of newton iterates (xn(:,end) is the final answer)

xn = zeros(length(x0),Nmax);    % pre allocate space for iterates
xn(:,1) = x0 - J(x0)\f(x0);

for i = 2:Nmax
    xn(:,i) = xn(:,i-1) - J(xn(:,i-1))\f(xn(:,i-1));
    
    if norm(f(xn(:,i))) < tau
        break
    end
end
xn = xn(:,1:i);
end
