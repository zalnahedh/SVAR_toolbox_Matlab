function output = structToIrfRecursive(x,info)

% --Description---------------------------------------------------------------
% structToIrfRecursive is mapping from vectorized structural representation   % x = [vec(A0); vec(Aplus)] to vectorized IRF representation                   % output=[vec(L0);vec(L1);...;vec(Lp);vec(C)] using recursive        % formulation of impulse response functions 
%---------------------------------------------------------------------------        

%  SVAR model:
%  y(t)'*A(0) = z(t)'*C + y(t-1)'*A(1) + ... + y(t-p)'*A(p) + epsilon(t)
%
%  y(t) -         (n x 1) endogenous variables
%  epsilon(t) -   (n x 1) exogenous shocks
%  z(t) -         (k x 1) exogenous or deterministic variables
%  A(i) -         (n x n) matrix
%  C -            (k x n) matrix 
%
%  structural representation  - A0, Aplus
%
%    A0 = A(0)                      [matrix (n x n)]
%    Aplus = [A(1); ... A(p); C]    [matrix (m x n) where m=n*p+k] 
%    x = [vec(A0); vec(Aplus)]      [vector (nx x 1) where nx=n(n+m)]
%
% 
%  IRF representation – L0, L+=[L1;L2;...;Lp;C]
%
%    L0 = inv(A(0))'
%    Lk = [L(k-1)*A(1)' + ... + L0*A(k)']*inv(A(0)')  1 <= k <= p
%
%    output = [vec(L0);vec(L1);...;vec(Lp);vec(C)]
%
%----------------------------------------------------------------------------
 
n=info.nvar;
p=info.nlag;
 
n_squared=n*n;
nx=numel(x); 
m=nx/n-n;
k=m-n*p;

%extract matrix A0 from vectorized implementation x
A0=reshape(x(1:n_squared),n,n);
%extract matrix Aplus from vectorized implementation x
Aplus=reshape(x(n_squared+1:end),m,n);
 
%extract matrices A(1),...,A(p) from matrix Aplus
A=cell(p,1);
for i=1:p
	%extract submatrix A{i} from matrix Aplus
    	A{i}=Aplus((i-1)*n+1:i*n,:);
end

%extract matrix C from matrix Aplus
C=Aplus(p*n+1:end,:);
 
%IRF representation
L=cell(p,1);
 
%recursive calculations of IRFs
L0=inv(A0)';
output=[vec(L0)]; %vectorized implementation
for k=1:p
    X=L0*A{k}';
    for j=1:(k-1)
        X=X+L{k-j}*A{j}';
    end
    L{k}=X*L0;
    output=[output;vec(L{k})];
end
output=[output;vec(C)];
end
