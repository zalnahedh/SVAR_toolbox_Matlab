function output = f_h_map_inv(x,info)

% --Description------------------------------------------------------------
% f_h_map_inv is mapping from vectorized orthogonal reduced-form      
% representation to vectorized structural representation (characterized by     
% transformation h)
% -------------------------------------------------------------------------
 
 
nvar=info.nvar;
m=info.m;
 
% Create matrices B, Sigma and Q from vectorized implementation x:
% Matrix B dimension is: (m x nvar)
% Matrix Sigma dimension is: (nvar x nvar)
% Matrix Q dimension is: (nvar x nvar)
B=reshape(x(1:m*nvar),m,nvar);
Sigma=reshape(x(m*nvar+1:m*nvar+nvar^2),nvar,nvar); 
Q=reshape(x(m*nvar+nvar^2+1:end),nvar,nvar);
 
% Calculation of structural representation: A0 and Aplus
A0=h_tilda(Sigma,info)\Q;
Aplus=B*A0;
 
% Vectorized implementation of structural representation:
output=[vec(A0); vec(Aplus)];
end

