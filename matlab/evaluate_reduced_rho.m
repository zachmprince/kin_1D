function rrho = evaluate_reduced_rho(phi,curr_time)
% Evaluate the reduced reactivity at the current time

global dat npar 

D    = assemble_stiffness(dat.cdiff   ,curr_time);
A    = assemble_mass(     dat.siga    ,curr_time);
NFIp = assemble_mass(     dat.nusigf_p,curr_time) / npar.keff;

tmp=NFIp-(D+A);

rho=npar.phi_adj' * tmp * phi;
rrho=rho/npar.K0;

end

