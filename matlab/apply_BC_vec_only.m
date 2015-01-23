function rhs=apply_BC_vec_only(rhs)

global dat npar

n=npar.n;

% apply BC (Dirichlet 0 for the mass matrix where Dirichlet is required!!!)
Dirichlet_nodes=[];
Dirichlet_val=[];

if(dat.bc.left.type==2)
    Dirichlet_nodes=[Dirichlet_nodes 1];
    Dirichlet_val=[Dirichlet_val dat.bc.left.C];
    if length(dat.bc.left.C) == 2
        Dirichlet_nodes=[Dirichlet_nodes (n+1)];
    end
end
if(dat.bc.rite.type==2)
    Dirichlet_nodes=[Dirichlet_nodes n];
    Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
    if length(dat.bc.left.C) == 2
        Dirichlet_nodes=[Dirichlet_nodes (n+n)];
    end
end

% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
    rhs(id)=bcval;         % put the constrained value in the rhs
end

return
end