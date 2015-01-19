function [F]=fem5_ss_residual(curr_time,T)
% assemble the residual and apply BC

% make the problem-data a global variable
global dat npar

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=npar.ndofs;
% allocate memory
F=zeros(n,1);

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
[xq,wq] = GLNodeWt(2*(porder+1));
% initialize local residual vector
local_res=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

% definition of the weak form:
% int_domain (grad u D grad b) - int_bd_domain (b D grad u n) ...
%      + int_domain( XSa u b) = int_domain (b rhs)

% loop over elements
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % jacobian of the transformation to the ref. element
    Jac=(x1-x0)/2;
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
    % T is of length porder+1, 2/dx is the 1d jacobian
    local_T    = b(:,:) * T(gn(iel,:));
    local_dTdx = dbdx(:,:) * T(gn(iel,:));
    d=dat.diff(curr_time,x,npar.scale*local_T)/Jac;
    s=dat.siga(curr_time,x,npar.scale*local_T)*Jac;
    q=dat.esrc(curr_time,x,npar.scale*local_T)*Jac/npar.scale;
    % compute local residual
    for i=1:porder+1
        local_res(i) =  dot(s.*wq.*b(:,i)    , local_T(:)) ...
            + dot(d.*wq.*dbdx(:,i) , local_dTdx(:)) ...
            - dot(q.*wq, b(:,i));
    end
    % assemble
    F(gn(iel,:)) = F(gn(iel,:)) + local_res;
end

% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
switch dat.bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        F(1)=F(1)-dat.bc.left.C;
    case 1 % Robin
        F(1)=F(1)+1/2*T(1);
        F(1)=F(1)-2*dat.bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val dat.bc.left.C];
end
switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        F(n)=F(n)-dat.bc.rite.C;
    case 1 % Robin
        F(n)=F(n)+1/2*T(n);
        F(n)=F(n)-2*dat.bc.rite.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes n];
        Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
end
% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
    F(id)=T(id)-bcval;         % put the constrained value in the rhs
end

% % % apply Mass matrix here, for now
% % F=dat.Mass\F;

F=-F; % !!!!!!!!!!!!!!!

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

