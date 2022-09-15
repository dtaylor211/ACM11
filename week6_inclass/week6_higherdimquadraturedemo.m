% ACM 11 Demo Script
% Numerical integration in higher dimensions

% define three functions in 1D, 2D, and 3D to integrate numerically
f1 = @(x) sin(x.^2);
f2 = @(x,y) sin(x.^2 + y.^2);
f3 = @(x,y,z) sin(x.^2 + y.^2 + z.^2);

% use MATLAB numerical integrators to calculate reference integration value
exact1 = integral(f1,-1,1);
exact2 = integral2(f2,-1,1,-1,1);
exact3 = integral3(f3,-1,1,-1,1,-1,1);


% use tensor gauss-legendre quadrature rules to integrate the three
% functions on the [-1,1] interval/square/cube
Nmax = 9;
err = zeros(Nmax,3);
for N = 1:Nmax
    
    % one-dimensional quadrature
    [points,weights] = qrule(N);
    estimate1 = weights*f1(points');
	err(N,1) = abs(exact1 - estimate1);
    
    % two-dimensional quadrature
    [xx,yy] = meshgrid(points,points);
    ww = weights.*weights';
    ff = f2(xx,yy);
    estimate2 = sum(sum(ww.*ff));
    err(N,2) = abs(exact2-estimate2);
    
    % three-dimensional quadrature
    [xx,yy,zz] = ndgrid(points,points,points);
    [w,v,u] = ndgrid(weights,weights,weights);
    estimate3 = sum(sum(sum( w.*v.*u.*f3(xx,yy,zz))));
    err(N,3) = abs(exact3-estimate3);
end

% plot errors vs. number of quadrature points
figure(1); clf
loglog(1:Nmax,err(:,1)); hold on
loglog((1:Nmax).^2, err(:,2));
loglog((1:Nmax).^3, err(:,3));
xlabel('# function evaluations','FontSize',14)
ylabel('Numerical integration error','FontSize',14)
legend('1D integral','2D integral','3D integral','FontSize',14)

% The point of this demo is to notice that to reach a given error level,
% the number of quadrature points required increases exponentially with the
% dimension. E.g. you get similar errors with 10 points for 1D functions,
% 100 points for 2D functions, and 1000 points for 3D functions. Many
% advanced computational methods require the estimation of integrals in
% higher dimensions, and traditional quadrature becomes very expensive in
% those settings.