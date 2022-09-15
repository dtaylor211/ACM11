function [opt_path, conv_path] = Rosenbrock_newton_bisec(x0, y0, alpha_UB, epsilon, Nmax)
    xn = x0;
    yn = y0;
    opt_path = zeros(2, Nmax+1);
    conv_path = zeros(1, Nmax+1);
    for n = 1 : Nmax
        % call rosenbrock.m at (xn,yn) to get the current fn, gn, and Hn
        [fn,gn,Hn] = rosenbrock(xn,yn);
        
        % store fn in nth entry of conv_path, and store xn,yn in nth column
        % of opt path
        conv_path(n) = fn;
        opt_path(:,n) = [xn,yn];
        
        % check to see if norm of gradient is below epsilon tolerance - if
        % so, break the for loop.
        if norm(gn) < epsilon
            break
        end
        
        
        %% use bisection algorithm to find step size
        % define \hat{g}, the Newton update direction, as defined in the
        % main file - note that this is NOT the same as \hat{g} for
        % gradient descent
        g_hat = (Hn\gn) / norm(Hn\gn);
        
        % define q(alpha) for current iteration as described in main file
        q = @(alpha) -[-2*(1-(xn-alpha*g_hat(1)))-4*(xn-alpha*(g_hat(1)))*(yn-alpha*g_hat(2) ...
            -(xn-alpha*g_hat(1))^2), 2*((yn-alpha*g_hat(2))-(xn-alpha*g_hat(1)).^2)]...
            *inv([2-4*(yn-alpha*g_hat(2))+12*(xn-alpha*g_hat(1)).^2, -4*(xn-alpha*g_hat(1)); -4*(xn-alpha*g_hat(1)), 2]) * g_hat;
       
        % call bisection routine on q(alpha) to find optimal step size
        [iterates,] = bisection(q,0,alpha_UB,10^(-10),40);
        alpha_o = iterates(end);

        % take Newton step with optimal step size alpha (you must use the
        % normalized Newton direction here)
        xn = xn - alpha_o * g_hat(1); % replace with update
        yn = yn - alpha_o * g_hat(2); % replace with update

    end
    opt_path(:, n+1) = [xn; yn];
    conv_path(:, n+1) = rosenbrock(xn,yn);
    opt_path = opt_path(:, 1:n+1);
    conv_path = conv_path(:, 1:n+1);
end