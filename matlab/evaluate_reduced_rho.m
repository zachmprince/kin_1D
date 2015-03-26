function rho = evaluate_reduced_rho(phi,curr_time)
% Evaluate the reduced reactivity at the current time

global dat npar 

D    = assemble_stiffness(dat.cdiff   ,curr_time);
A    = assemble_mass(     dat.siga    ,curr_time);
NFI = assemble_mass(     dat.nusigf,curr_time) / npar.keff;

tmp=NFI-(D+A);

rho=(npar.phi_adj' * tmp * phi) / npar.K0;

end