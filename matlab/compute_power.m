function pow = compute_power(fct_ptr,curr_time,phi)

% make the problem-data a global variable
global npar

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=npar.n;
% allocate memory
vec=zeros(n,1);
% initialize local matrices
v=zeros(porder+1,1);
% shape set
b=npar.b;
% quadrature stuff
xq=npar.xq;
wq=npar.wq;


pow=0;
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
    % 2/dx is the 1d jacobian
    % old: s=fct_ptr(curr_time,x)*Jac;
    imat=npar.elem_to_mat(iel);
    s=evaluate_material_prop(fct_ptr{imat},curr_time,x)*Jac/npar.keff;

    % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
    % local flux
    local_phi = b(:,:) * phi(gn(iel,:));
    % integrate over cell
        pow = pow + dot(s.*wq ,local_phi);
end

return
end

