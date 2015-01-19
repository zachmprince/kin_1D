function M=assemble_transient_operator(curr_time)

global dat npar

n   = npar.n;
nnz = npar.nnz;
M = sparse(2*n,2*n,2*(nnz+n));

D    = assemble_stiffness(dat.diff,curr_time);
A    = assemble_mass(dat.siga,curr_time);
NFIp = assemble_mass(dat.nusigf_prompt,curr_time);
NFId = assemble_mass(dat.nusigf_delayed,curr_time);

tmp=NFIp-(D+A); 
tmp=apply_BC_mat_only(tmp,true);
NFId=apply_BC_mat_only(NFId,true);

L=dat.lambda*speye(n);
Ldiag=apply_BC_mat_only(L,true);
Loffd=apply_BC_mat_only(L,true);



M(1:n    ,1:n    ) = tmp;
M(1:n    ,n+1:2*n) =  Loffd;
M(n+1:2*n,1:n    ) = NFId;
M(n+1:2*n,n+1:2*n) = -Ldiag;

return
end