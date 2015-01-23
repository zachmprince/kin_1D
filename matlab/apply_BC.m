function [A,rhs]=apply_BC(A,rhs,zero_out)

global dat npar

n=npar.n;

% apply BC (Dirichlet 0 for the mass matrix where Dirichlet is required!!!)
Dirichlet_nodes=[];
Dirichlet_val=[];

if(dat.bc.left.type==2)
    Dirichlet_nodes=[Dirichlet_nodes 1];
    Dirichlet_val=[Dirichlet_val dat.bc.left.C];
end
if(dat.bc.rite.type==2)
    Dirichlet_nodes=[Dirichlet_nodes n];
    Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
end

if zero_out
    diag_value=0;
else
    diag_value=1;
end

% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
    rhs(1:n)=rhs(1:n)-bcval*A(1:n,id);  % modify the rhs using constrained value
    A(id,:)=0; % set all the id-th row to zero
    A(1:n,id)=0; % set all the id-th column to zero (symmetrize A)
    A(id,id)=diag_value;            % set the id-th diagonal to unity
    rhs(id)=bcval;         % put the constrained value in the rhs
end

return
end