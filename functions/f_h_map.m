function output = f_h_map(x,info)

% --Description---------------------------------------------------------------
% f_h_map is mapping from vectorized structural representation to vectorized   % orthogonal reduced-form representation (characterized by transformation h)
% ----------------------------------------------------------------------------
 
nvar=info.nvar;
m=info.m;
 
% Create matrices A0 and Aplus from vectorized implementation x:
% Matrix A0 dimension is: (nvar x nvar)
% Matrix Aplus dimension is: (m x nvar)
A0=reshape(x(1:nvar*nvar),nvar,nvar);
Aplus=reshape(x(nvar*nvar+1:end),m,nvar);
 
% Calculation of orthogonal reduced-form: B, Sigma and Q
B=Aplus/A0;
Sigma=inv(A0*A0');
Q=h_tilda(Sigma,info)*A0;
 
% Vectorized implementation of orthogonal reduced-form representation:
output=[vec(B); vec(Sigma); vec(Q)];
end

