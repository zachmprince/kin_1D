function kinetics_driver

% clear all;
close all; clc;
global dat npar

% % get current directory
% s=what; cur_dir=s.path; addpath(genpath(cur_dir))

% select problem
pbID=1;
problem_init(pbID);

% compute eigenmode
curr_time=0;
[phi,keff]=steady_state_eigenproblem(curr_time);
plot(npar.x_dofs,phi)

% hold all
% [phi1,keff1]=steady_state_eigenproblem(1);
% plot(dat.x_dofs,phi1)
%
% [phi2,keff2]=steady_state_eigenproblem(3);
% plot(dat.x_dofs,phi2)
%
% [keff keff1 keff2]'

C=kinetics_init(phi,curr_time);
% npar.K0
% npar.Pnorm
u=[phi;C]; u0=u;

dt=0.01;
Pnorm=npar.Pnorm;
for it=1:500
    time_end=it*dt;
    fprintf('time end = %g \n',time_end);

    TR=assemble_transient_operator(time_end);
    M =assemble_time_dependent_operator(time_end);
    % M(unew-uold)/dt=TR.unew
    rhs=M*u;
    u = (M-dt*TR)\rhs;
    plot(dat.x_dofs,u(1:npar.n));drawnow

    POW=assemble_load(dat.nusigf,time_end);
    Pnorm(it+1)=dot(POW,u(1:npar.n));


end

figure(2);
plot(Pnorm/Pnorm(1),'+-')


% rm path that was added
rmpath(dat.full_dir);

return
end