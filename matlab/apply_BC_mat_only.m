function A=apply_BC_mat_only(A,zero_out)

global dat npar

n=npar.n;

% apply BC (Dirichlet 0 for the mass matrix where Dirichlet is required!!!)
Dirichlet_nodes=[];
if(dat.bc.left.type==2)
    Dirichlet_nodes=[Dirichlet_nodes 1];
end
if(dat.bc.rite.type==2)
    Dirichlet_nodes=[Dirichlet_nodes n];
end

if zero_out
    diag_value=0;
else
    diag_value=1;
end

% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    A(id,:)=0; % set all the id-th row to zero
    A(:,id)=0; % set all the id-th column to zero (symmetrize A)
    A(id,id)=diag_value;            % set the id-th diagonal to unity or 0
end

return
end