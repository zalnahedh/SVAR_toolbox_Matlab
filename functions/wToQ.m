function output = wToQ(x,info)
%----------------------------------------------------------------------
% input:  x=[vec(B);vec(Sigma);vec(w1,..,wn)]
% output: orthogonal Q matrix
%----------------------------------------------------------------------

nvar=info.nvar;
m=info.m;
W_matrix=info.W_matrix;

% extract reduced form paramteres from input vector x
extract_reduced_param=x(1:(m+nvar)*nvar);
% extract vectors w from input vector x ( w1,w2,...,wn)
w=x((m+nvar)*nvar+1:end);

% Transform (Beta,Sigma,In) to structural representation z
z=f_h_map_inv([extract_reduced_param;vec(eye(nvar))],info);

% Evaluate ZF for given structural representation z:
ZF=info.ZF(z,info);

% Initialize the orthogonal matrix Q:
Q=zeros(nvar,nvar);

k=0;
for j=1:nvar
    s=size(W_matrix{j},1);
    wj=w(k+1:k+s);
    Mj_tilde=[Q(:,1:j-1)'; ZF{j}; W_matrix{j}];
    [K,R]=qr(Mj_tilde');    
    for i=nvar-s+1:nvar
        if (R(i,i) < 0) 
            K(:,i)=-K(:,i);
        end
    end 
    Kj=K(:,nvar-s+1:nvar);
    Q(:,j)=Kj*wj;
    k=k+s;
end
output=reshape(Q,nvar*nvar,1);
end
