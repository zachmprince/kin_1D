function [M,P]=assemble_steady_state_operator(curr_time)

global dat

D  = assemble_stiffness(dat.cdiff,curr_time);
A  = assemble_mass(dat.siga,curr_time);
M=D+A;

P  = assemble_mass(dat.nusigf,curr_time);

return
end