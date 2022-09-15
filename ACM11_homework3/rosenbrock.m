function [f,g,H] = rosenbrock(x,y)
% "nargout" is the number of output arguments requested
% if we call this function and ask for 1 output, it will just return f
% if we ask for 2 outputs, we will get f and g
% if we ask for 3 outputs, we get f, g, and H

f = (1-x).^2 + (y-x.^2).^2;

if nargout >= 2 
    g1 = -2*(1-x) - 4*x*(y-x^2);
    g2 = 2*(y-x^2);
    g = [g1; g2]; % fill in correct gradient expression here
end

if nargout == 3
    h11 = 2 - 4*y + 12*x^2;
    h12 = -4*x;
    h21 = -4*x;
    h22 = 2;
    H = [h11, h12; h21, h22]; % fill in correct Hessian expression here
end
