function [beta] = drawBetaCondSigma(info,beta0,omega0,sigma)

q=info.nvar*(info.nvar*info.nlag+info.nex+info.cte);

% first version:
% omega0Inv=invpd(omega0);
% sigmaInv=invpd(sigma);
% omega1=invpd(omega0Inv+kron(sigmaInv,(info.X)'*info.X));
% previous steps have numerical issues

% omega0 is diagonal matrix, inverse is easily calculated
omega0Inv=diag(1./diag(omega0));

% if A is SPD matrix then A= R' * R 
% chol(A) is R (upper trangular matrix)
% (R)^(-1) = R \ I  => (A)^(-1)=(R' * R )^(-1)=[(R)^(-1)] * [(R)^(-1)]'
% If we are not sure that matrix is SPD, first find nearest SPD matrix
% Calculate Cholesky decomposition of SPD matrix to calculate inverse

help=chol(findNearestSpd(sigma));
invhelp=help\speye(info.nvar);
sigmaInv=invhelp*invhelp';

omega1Inv=omega0Inv+kron(sigmaInv,(info.X)'*info.X);
help=chol(findNearestSpd(omega1Inv));
invhelp=help\speye(q);
omega1=invhelp*invhelp';

beta1=omega1*(omega0Inv*beta0+kron(sigmaInv,(info.X)')*vec(info.Y)); 
omega1=(omega1+omega1')/2; %symmetric
beta=beta1+(randn(1,q)*chol(omega1))';

end

