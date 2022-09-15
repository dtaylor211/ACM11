function xn = bisection(f,xa,xb,tau,Nmax)
% Solves f(x) = 0 via bisection
% Inputs: 
%    f     function to find root of
%    xa,xb initial search interval bounds
%    tau   tolerance for convergence
%    Nmax  max number of iterations
% Outputs:
%    xn    sequence of iterates

xn = zeros(length(xa),Nmax);        % pre allocate space for iterates

assert(sign(f(xa))~=sign(f(xb)))
for i = 1:Nmax
    xn(:,i) = 0.5*(xa+xb);
    if norm(f(xn(:,i))) < tau
        break
    elseif sign(f(xa))==sign(f(xn(:,i)))
        xa = xn(:,i);
    else
        xb = xn(:,i);
    end
end
xn = xn(:,1:i);
end