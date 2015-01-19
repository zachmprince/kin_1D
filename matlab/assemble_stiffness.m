function A = assemble_stiffness(fct_ptr,curr_time)

% make the problem-data a global variable
global npar

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=npar.n;
% ideally, we would analyze the connectivity to determine nnz
nnz=npar.nnz;
% allocate memory
A=spalloc(n,n,nnz);
% initialize local matrices
a=zeros(porder+1,porder+1);
% shape set
dbdx=npar.dbdx;
% quadrature stuff
xq=npar.xq;
wq=npar.wq;

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
    d=fct_ptr(curr_time,x)/Jac;

    % assemble
    for i=1:porder+1
        for j=1:porder+1
            a(i,j)= dot(d.*wq.*dbdx(:,i)    , dbdx(:,j));
        end
    end
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + a;
end

return
end

