function uT = forwardEulerLinear(u0,kappa,A,dt,T)
    N = length(u0);
    K = ceil(T/dt);
    uT = u0;
    for j = 0:K
        uT = (eye(N) + kappa*dt*A)*uT;
    end
end

