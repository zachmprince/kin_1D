function [rho_MGT,beff_MGT,prec]=compute_prke_parameters(curr_time,shape,C)

global dat npar

phi_adjoint = npar.phi_adj;

D    = assemble_stiffness(dat.diff,curr_time);
A    = assemble_mass(dat.siga,curr_time);
NFI  = assemble_mass(dat.nusigf,curr_time);
NFId = assemble_mass(dat.nusigf_delayed,curr_time);
IV   = assemble_mass(dat.inv_vel,curr_time);

rho  = phi_adjoint' * (NFI-D_A) * shape;
beff = phi_adjoint' * NFId * shape;
MGT  = phi_adjoint' * IV * shape;

rho_MGT  = rho/MGT;
beff_MGT = beff/MGT;

L=dat.lambda*speye(n);

prec = phi_adjoint' * L * C;


return
end