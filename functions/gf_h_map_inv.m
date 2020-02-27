function y =gf_h_map_inv(x,info)
%----------------------------------------------------------------------
% input:  x=[vec(B);vec(Sigma);vec(w1,..,wn)]
% output: [vec(A0);vec(Aplus)] 
% Output is structural representation that satisfies zero restrictions
%----------------------------------------------------------------------
% g_map_inv            : (Beta,Sigma,w1,...,wn) -> (Beta,Sigma,q1,...,qn)
% f_h_map_inv          : (Beta,Sigma,q1,...,qn) -> (A0,Aplus)
% f_h_map_inv o g_map_inv  : (Beta,Sigma,w1,...,wn) -> (A0,Aplus) 
% f_h_map_inv o g_map_inv <=> (g_map o f_h_map)_inv  which we call 
% gf_h_map_inv 
%-----------------------------------------------------------------------
% gf_h_map_inv is mapping from (Beta,Sigma,W) to structural representation 
% (A0,Aplus) and by construction this structural representation satisfies 
% zero restrictions
%-----------------------------------------------------------------------
% It is important to note that the zero rerstrictions in the orthogonal 
% reduced-form parameterization are just linear restrictions on each column
% of the orthogonal matrix Q, conditional on the reduced form parameters!
% Zj*F(f_h_inv(Beta,Sigma,Q))ej = Zj* F(f_h_inv(Beta,Sigma,In))*Q*ej=0
%------------------------------------------------------------------------

nvar=info.nvar;
m=info.m;
W_matrix=info.W_matrix;

% extract reduced form paramteres from input vector x
extract_reduced_param=x(1:(m+nvar)*nvar);
% extract vectors w from input vector x ( w1,w2,...,wn)
w=x((m+nvar)*nvar+1:end);

% Transform (Beta,Sigma,In) to structural representation z
z=f_h_map_inv([extract_reduced_param;vec(eye(nvar))],info);

%structural representation in notation (A0,Aplus):
A0=reshape(z(1:nvar*nvar),nvar,nvar);
Aplus=reshape(z(nvar*nvar+1:end),m,nvar);

% Evaluate ZF for given structural representation z:
ZF=info.ZF(z,info);

% Initialize the orthogonal matrix Q:
Q=zeros(nvar,nvar);

k=0;
for j=1:nvar
    z=size(W_matrix{j},1);
    wj=w(k+1:k+z);
    % Mj_tilde=[Mj' wj']'
    % where: Mj=[Qj-1 ZjF(f_h_map_inv(Beta,Sigma,In)))' ]'
    Mj_tilde=[Q(:,1:j-1)'; ZF{j}; W_matrix{j}];
    [K,R]=qr(Mj_tilde');    
    for i=nvar-z+1:nvar
        if (R(i,i) < 0) 
            K(:,i)=-K(:,i);
        end
    end 
    Kj=K(:,nvar-z+1:nvar);
    %calculate j-th column of orthogonal matrix Q:
    Q(:,j)=Kj*wj;
    k=k+z;
end

% final output is structural representation that satisfies zero
% restrictions:
y=[vec(A0*Q); vec(Aplus*Q)];

end

