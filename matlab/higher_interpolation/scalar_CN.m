clear all; close all; clc;

% we solve du/dt = a(t).u(t) + q(t) with u(t=0)=u0
% a(t) is an arbitrary function of time
% when employing a user-chosen exact solution, u=ue, then this
% automatically sets q(t)= u'-a.u. otherwise, q(t)=0

% end of simulation time
tend=1;
% coefficients for time dependence of exact solution, if selected
coef = [1 1 1 1]; % coef for t^0, t^1, t^2, t^3
a_const=2;
% a(t) choice
a_type = 3;
% solution type
solution_type = 3;

% tolerances for odesolvers
rtol = 1e-13; abso = 1e-13;
atol  = abso*ones(length([1]),1);
options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);

% pick function for a(t)
switch a_type
    case 0 % constant
        a = @(t)(0.*t+a_const);
        dadt = @(t)(0.*t);
    case 1 % linear
        a = @(t)(t);
        dadt = @(t)(0.*t+1.);
    case 2 % decaying exponential
        a = @(t)(exp(-t));
        dadt = @(t)(-exp(-t));
    case 3 % increasing exponential
        a = @(t)(exp(t));
        dadt = @(t)(exp(-t));
    otherwise
        error('unknown type for function a(t)');
end


% pick function for q(t)
switch solution_type
    case 0 % no manufactured solution
        % set q(t)
        q = @(t) (0.*t);
        % impose initial value
        u0 = 1;
        % solve the ODE to get the final time reference value
        [tt,uu]=ode15s(@scalar_ssres,[0 tend],u0,options,a,q);
        % extract reference value
        u_ref = uu(end);
        % plot
        plot(tt,uu); hold all;

    case 2 % u_exact is quadratic
        exact = @(t)(coef(1)+coef(2)*t+coef(3)*t.^2);
        dexactdt = @(t)(coef(2)+2*coef(3)*t);
        % set q(t)
        q = @(t) (dexactdt(t) - a(t)*exact(t));
        % reference value
        u_ref = exact(tend);
        % plot
        tt=linspace(0,tend);
        plot(tt,exact(tt)); hold all;
        % impose initial value
        u0 = exact(0);
        % solve the ODE to get the final time reference value
        [tt,uu]=ode15s(@scalar_ssres,[0 tend],u0,options,a,q);
        % plot
        plot(tt,uu);

    case 3 % u_exact is cubic
        exact = @(t)(coef(1)+coef(2)*t+coef(3)*t.^2+coef(4)*t.^3);
        dexactdt = @(t)(coef(2)+2*coef(3)*t+3*coef(4)*t.^2);
        % set q(t)
        q = @(t) (dexactdt(t) - a(t)*exact(t));
        % reference value
        u_ref = exact(tend);
        % plot
        tt=linspace(0,tend);
        plot(tt,exact(tt)); hold all;
        % impose initial value
        u0 = exact(0);
        % solve the ODE to get the final time reference value
        [tt,uu]=ode15s(@scalar_ssres,[0 tend],u0,options,a,q);
        % plot
        plot(tt,uu);

    otherwise
        error('unknown type for solution');
end


%%% convergence study
n_runs=10;
n_steps = 2.^linspace(1,n_runs,n_runs);

for i_run=1:n_runs
    
    nbr_steps = n_steps(i_run);
    dt = tend/nbr_steps;
    
    sol=zeros(nbr_steps+1,1);
    sol(1)=u0;
    
    for it=1:nbr_steps
        time0 = (it-1)*dt;
        time1 = time0 + dt;
        % CN: (unew - uold)/dt = 0.5 ( a(old)u(old)+q(old) + a(new)u(new)+q(new) )
        deno = 1-dt/2*a(time1);
        sol(it+1) = ( sol(it) + dt/2*(scalar_ssres(time0,sol(it),a,q) + q(time1)) )/deno;
    end
    plot(linspace(0,tend,nbr_steps+1),sol);
    conv_sol(i_run) = abs(sol(end)-u_ref);
    conv_dt(i_run) = dt;
    
end

figure(2)
loglog(conv_dt,conv_sol,'+-'); hold all;
% error = C dt^2. Compute C for the smallest dt value
C = conv_sol(end)/conv_dt(end)^2;
loglog(conv_dt,C*conv_dt.^2);
    

[conv_dt' conv_sol']