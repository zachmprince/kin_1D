function [A]=compute_ss_precond(curr_time,T)
% assemble the P matrix, with various options:
% optP=0: Identity
% optP=1: compute k(T) at quad points
% optP=1a: compute k(T) at quad points using T^(0)
% optP=2: compute kave per element
% optP=2a: compute kave per element using T^(0)
% the options with a are actually managed outside of this function

% make the problem-data a global variable
global dat npar

optP=npar.prec_opt;

if(optP==0)
    A=speye(length(T));
    return
end

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+1)*nel; %this is an upperbound, not exact
% n: linear system size
n=npar.ndofs;
% allocate memory
A=spalloc(n,n,nnz);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
k=m;
local_mat=m;

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
if(optP==1)
    poly_max=2*(porder+1);
    [xq,wq] = GLNodeWt(poly_max);
    % store shapeset
    [b,dbdx] =feshpln(xq,porder);
else
    poly_max=porder+1;
    [xq,wq] = GLNodeWt(poly_max);
    % store shapeset
    [b,dbdx] =feshpln(xq,porder);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            m(i,j)= dot(wq.*b(:,i)    , b(:,j));
            k(i,j)= dot(wq.*dbdx(:,i) , dbdx(:,j));
        end
    end
    % increase quadrature order to compute kave over each element
    poly_max=2*(porder+1);
    [xq,wq] = GLNodeWt(poly_max);
    % store shapeset
    [b,dbdx] =feshpln(xq,porder);
end

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
    if(optP==1)
        % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
        local_T    = b(:,:) * T(gn(iel,:));
        d=dat.diff(curr_time,x,npar.scale*local_T);
        s=dat.siga(curr_time,x,npar.scale*local_T);
        % compute local matrices
        for i=1:porder+1
            for j=1:porder+1
                m(i,j)= dot(s.*wq.*b(:,i)    , b(:,j));
                k(i,j)= dot(d.*wq.*dbdx(:,i) , dbdx(:,j));
            end
        end
        local_mat = m*Jac + k/Jac;
    else
        % first compute ave_cond
        % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
        local_T    = b(:,:) * T(gn(iel,:));
        d=dat.diff(curr_time,x,npar.scale*local_T);
        ave_cond = dot(wq, d)/sum(wq);
        % then compute local mat
        local_mat =  0*m*Jac + ave_cond*k/Jac;
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + local_mat;
end

% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
switch dat.bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
%         rhs(1)=rhs(1)+dat.bc.left.C;
    case 1 % Robin
        A(1,1)=A(1,1)+1/2;
%         rhs(1)=rhs(1)+2*dat.bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val dat.bc.left.C];
end
switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
%         rhs(n)=rhs(n)+dat.bc.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n)+1/2;
%         rhs(n)=rhs(n)+2*dat.bc.rite.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes n];
        Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
end
% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
%     rhs=rhs-bcval*A(:,id);  % modify the rhs using constrained value
    A(id,:)=0; % set all the id-th row to zero
    A(:,id)=0; % set all the id-th column to zero (symmetrize A)
    A(id,id)=1;            % set the id-th diagonal to unity
%     rhs(id)=bcval;         % put the constrained value in the rhs
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
