function output = irfToStructRecursive(x,info)

% --Description---------------------------------------------------------------
% irfToStructRecursive is mapping from vectorized IRF representation          % x=[vec(L0);vec(L1);...;vec(Lp);vec(C)] to vectorized structural    % representation output = [vec(A0); vec(Aplus)] using recursive                  % formulation of impulse response functions 
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
%  IRF representation – L0, L+=[L1;L2;...;Lp;C]%
%
%  structural representation  - A0, Aplus = [A(1); ... A(p); C
%
%----------------------------------------------------------------------------
n=info.nvar;
p=info.nlag;
 
n_squared=n*n;
nx=numel(x); 
m=nx/n-n;
k=m-n*p;
 
%extract matrix L0 from vectorized implementation x
L0=reshape(x(1:n_squared),n,n);
%extract matrix Lplus from vectorized implementation x
Lplus=reshape(x(n_squared+1:end),m,n);
 
%extract matrices L1,...,Lp from matrix Lplus
L=cell(p,1);
for i=1:p
    %extract submatrix Li from matrix Lplus
        L{i}=Lplus((i-1)*n+1:i*n,:);
end
 
%extract matrix C from matrix Lplus
C=Lplus(p*n+1:end,:);
 
%structural representation
A=cell(p,1);
 
%recursive calculations of structural representation:
A0=inv(L0)';
output=[vec(A0)]; %vectorized implementation
for k=1:p
    X=A0*L{k}';
    for j=1:(k-1)
        X=X-A{k-j}*L{j}';
    end
    A{k}=X*A0;
    output=[output;vec(A{k})];
end
output=[output;vec(C)];
end
