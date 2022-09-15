n = 11;
x = linspace(-1,1,n)';
y = 1./(1+25*x.^2); 

fspline = spline(x,y); % Define spline interpolant function

xplot = linspace(-1,1,1000); % Define plotting domain
yexact = 1./(1+25*xplot.^2); 
ysplineinterp = ppval(fspline,xplot);
err_interp = max(abs(yexact-ysplineinterp));
disp(err_interp)

figure(5); clf
plot(xplot,ysplineinterp,'linewidth',2); hold on
plot(xplot,yexact,':'); 
plot(x,y,'kx','MarkerSize',10,'linewidth',2); hold on;
for i = 1:length(x)
   plot([x(i),x(i)],[0 1],'Color',[0 0 0 0.2],'LineWidth',1); 
end

hold off
xlabel('$x$','interpreter','LaTex'); 
ylabel('$y$ ','interpreter','LaTex'); 
legend('Spline interpolant','Exact function','Data points','Location','nw')
set(gca,'fontsize',10); 
title('Spline interpolation P(x) and exact function f(x)')

figure(6); clf
plot(xplot,abs(yexact-ysplineinterp),'linewidth',2); hold on
plot(x,zeros(size(x)),'kx','MarkerSize',10); hold off
xlabel('$x$','interpreter','LaTex'); 
ylabel('Absolute Error ','interpreter','LaTex'); 
set(gca,'fontsize',10); 
title('Absolute difference between P(x) and f(x)')