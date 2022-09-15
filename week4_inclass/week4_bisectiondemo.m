% ACM 11 Week 4: Bisection Demo
% 
% This file sets up and calls a bisection root finder to solve f(x) = 0 for
% two examples:
%   Ex 1: f(x) = cos(x) - x
%   Ex 2: f(x) = sin(x) - x^3
% Make sure the example you want to run is uncommented, and the example you
% don't want to run is commented out.
%
% Two figures are generated: in Figure 101, an animation of the bisection
% search is displayed. In Figure 1, we plot the guess at each step of the
% iteration and the distance from the final solution vs. iteration #.

clear
close all

tau = 1e-6; % tolerance for convergence
Nmax = 100; % max # iterations

% Example 1
f = @(x) cos(x) - x;
xl = -pi;
xr = pi;
[xn,intervals] = bisection(f,xl,xr,tau,Nmax);

% Example 2
% f = @(x) sin(x) - x.^3;
% xl = -2;
% xr = -0.5;
% [xn,intervals] = bisection(f,xl,xr,tau,Nmax);

figure(1); clf
subplot(2,1,1)
plot(xn,'k'); grid on; hold on; plot(xn,'ko')
title('Bisection method: iterates')
ylabel('x_n')
subplot(2,1,2)
semilogy(abs(xn-xn(end))); grid on
title('Distance away from final interate')
xlabel('Iteration number')

pseconds = 0.5; % # seconds to pause between displaying updates, can adjust
animate_bisect(f,xn,intervals,pseconds) 

function animate_bisect(f,xn,intervals,p)
bounds = max(max(abs(intervals)));
x = linspace(-bounds,bounds);

for i = 1:length(xn)
    figure(101); clf
    plot(x,f(x)); grid on; hold on
    xlabel(['x^{(n)}=',num2str(xn(i))],'FontSize',18)
    ylabel('f')
    plot(intervals(1,i),f(intervals(1,i)),'ro')
    plot(intervals(2,i),f(intervals(2,i)),'ro')
    plot(xn(i),f(xn(i)),'bo')
    diff = intervals(2,i) - intervals(1,i);
    title(['i = ',num2str(i),': interval size = ',num2str(diff)],'FontSize',18)
    plot([-bounds bounds],[0 0],'k','LineWidth',1)

    pause(p)    % pauses for p seconds
    drawnow;    % tells MATLAB to display the plot above now (usually MATLAB would run through the loop and just show us the last plot)
end
end