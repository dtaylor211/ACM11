% ACM 11 Week 4: Newton Demo
% 
% This file sets up and calls a Newton solver to solve f(x) = 0 for
% three examples:
%   Ex 1: f(x) = exp(x) - 4x^2 - 8;
%   Ex 2: f(x) = cos(x^2)
%   Ex 3: f(x) = (x-4)*(x+1)
%   Ex 4: f(x) = x^3 - 2x + 2
% Make sure the example you want to run is uncommented, and the examples 
% you don't want to run is commented out.
%
% Two figures are generated: in Figure 101, an animation of the Newton
% update is displayed. In Figure 1, we plot the guess at each step of the
% iteration and the distance from the final solution vs. iteration #.

clear; close all

% iteration parameters
Nmax = 20;
tau = 1e-6;

%% exponential function
f = @(x) exp(x) - 4*x.^2 -8;
J = @(x) exp(x) - 8*x;
x0 = 0;
xn = newton(f,J,x0,tau,Nmax);
fstr = 'f(x) = e^x - 4x^2 - 8';

%% cos(x^2)
% f = @(x) cos(x.^2);
% J = @(x) -sin(x.^2)*2*x;
% x0 = 0.2;   % can try 0.2, 1, 2 -- notice difference in final roots
% xn = newton(f,J,x0,tau,Nmax);        
% fstr = 'f(x) = cos(x^2)';

%% quadratic 
% f = @(x) (x-4).*(x+1);
% J = @(x) 2*x - 3;
% 
% x0 = 1.1; % can try -1, 1.1, 1.4, (1.5 is critical)
% xn = newton(f,J,x0,tau,Nmax);
% fstr = 'f(x) = x^2 - 3x - 4';

%% nonconvergent
% f = @(x) x.^3 - 2*x + 2;
% J = @(x) 3*x.^2 - 2;
% x0 = 0; % non convergence x0 = 0,1 -- but x0 = 0.5 would work
% xn = newton(f,J,x0,tau,Nmax);
% fstr = 'f = x^3 - 2x +2';

pseconds = 1;
animate_newton(f,J,xn,pseconds)

figure(1); 
subplot(2,1,1)
plot(xn); grid on
ylabel('Iterate')
title(['Newton Convergence, ',fstr])
subplot(2,1,2)
semilogy(abs(xn-xn(end))); grid on
xlabel('Iteration #')
ylabel('Error')

function animate_newton(f,J,xn,p)
fprintf('animating... ')
a = xn(end);
bounds = max(abs(xn-a));
xplt = linspace(a-bounds,a+bounds);
fplt = f(xplt);

ymin = min(min(fplt)-0.1*(max(fplt)-min(fplt)), -0.1);
ymax = max(fplt)+0.1*(max(fplt)-min(fplt));

for i = 1:length(xn)-1
    figure(101); clf
    plot(xplt,fplt); grid on; hold on
    xlabel(['x^{(n)}=',num2str(xn(i+1))],'FontSize',16)
    ylabel('f')
    
    plot(xn(i),f(xn(i)),'ko')
    tangent_line = J(xn(i))*(xplt-xn(i))+f(xn(i));
    plot(xplt,tangent_line,'k:')
    
    
    diff = (xn(i+1)-xn(end))./xn(end);
    title(['n = ',num2str(i),': rel. error = ',num2str(diff)],'FontSize',16)
    plot([a-bounds,a+bounds],[0 0],'k','LineWidth',1)
    ylim([ymin ymax])
    
    pause(0.5*p)
    plot([xn(i+1) xn(i+1)],[0,f(xn(i+1))],'r:')
    plot(xn(i+1),f(xn(i+1)),'ro')
    
    pause(0.5*p)
    drawnow;
    fprintf(' ... ')
end
fprintf('done.\n')
end