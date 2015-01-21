function [C]=kinetics_init(phi,curr_time)

% make the problem-data a global variable
global dat npar 

NFId = assemble_mass(dat.nusigf_delayed,curr_time);
NFId = apply_BC_mat_only(NFId,true);
C = NFId*phi/dat.lambda;

% verif
% M=assemble_transient_operator(curr_time);
% plot(M*[phi;C])

npar.phi_adj=phi;

IV = assemble_mass(dat.inv_vel,curr_time);
K0 = npar.phi_adj' * IV * phi;

POW=assemble_load(dat.nusigf,curr_time);
Pnorm=dot(POW,phi);

npar.K0 = K0;
npar.Pnorm = Pnorm;

return
end


