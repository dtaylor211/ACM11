% ACM 11 Week 4: Fixed Point Iteration Demo
% 
% This file sets up and calls a fixed point iteration to solve f(x) = 0 for
% two examples:
%   Ex 1: f(x) = cos(x) - x
%   Ex 2: f(x) = sin(x) - x^3
% Make sure the example you want to run is uncommented, and the example you
% don't want to run is commented out.
%
% Two figures are generated: in Figure 101, an animation of the bisection
% search is displayed. In Figure 1, we plot the guess at each step of the
% iteration and the distance from the final solution vs. iteration #.

clear; close all

% in this case the tolerance is in relative change from one guess to the
% next: if it is less than 0.1%, we say we have converged
tau = 1e-4; 

Nmax = 50;  % max iterations

% Example 1
x0 = 1;
f = @(x) cos(x) - x;
g = @(x) cos(x);

% Example 2
% x0 = pi/2;
% f = @(x) sin(x) - x.^3;
% g = @(x) sin(x).^(1/3);
% g = @(x) sin(x) -x.^3 + x;

xn = fpi(g,x0,tau,Nmax);

pseconds = 0.4;
animate_fpi(f,xn,pseconds)

figure(1); clf
subplot(2,1,1)
plot(real(xn),'k'); grid on; hold on
plot(real(xn),'ko')
xlabel('n')
ylabel('x_n')
title('FPI iterates')

subplot(2,1,2)
semilogy(abs(real(xn-xn(end))),'k'); grid on; hold on
semilogy(abs(real(xn-xn(end))),'ko')
xlabel('n')
ylabel('Distance from final answer')
title('FPI convergence')

function animate_fpi(f,xn,p)
    x = linspace(-pi,pi);
    for i = 1:length(xn)-1
        figure(101); clf
        plot(x,f(x)); grid on; hold on
        xlabel('x')
        ylabel('f')
        plot(xn(i),f(xn(i)),'ro')
        diff = (xn(i+1)-xn(i))./xn(i);
        title(['i = ',num2str(i),': \Delta_{rel} = ',num2str(diff)],'FontSize',18)
        pause(p)
        drawnow;
    end
end
