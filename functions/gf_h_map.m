function y =gf_h_map(x,info)
%----------------------------------------------------------------------------------
% input: [vec(A0);vec(Aplus)] structural representation 
% output:  x=[vec(B);vec(Sigma);vec(w1,..,wn)]
% If (A0,Aplus) satisfies the zero restrictions, then gf_h_map_inv(gf_h_map(A0,A+)) = (A0,A+)
%--------------------------------------------------------------------------------------
nvar=info.nvar;
m=info.m;
W_matrix=info.W_matrix;
dimSphere=info.dimSphere;

% Transform (A0,Aplus) to orthogonal reduced form representation z1:
z1=f_h_map(x,info);

% extract reduced form paramteres from vector z1
extract_reduced_param=z1(1:(m+nvar)*nvar);

% extract orthogonal Q matrix from orthogonal reduced form representation z1:
Q=reshape(z1(nvar*(m+nvar)+1:end),nvar,nvar);

% Transform (Beta,Sigma,In) to structural representation z2
z2=f_h_map_inv([extract_reduced_param;vec(eye(nvar))],info);

% Evaluate ZF for given structural representation z2:
ZF=info.ZF(z2,info);

% initialize vector w
w=zeros(dimSphere,1);
k=0;

for j=1:nvar
    z=size(W_matrix{j},1);
    Mj_tilde=[Q(:,1:j-1)'; ZF{j}; W_matrix{j}];
    [K,R]=qr(Mj_tilde');    
    for i=nvar-z+1:nvar
        if (R(i,i) < 0) 
            K(:,i)=-K(:,i);
        end
    end 
    Kj=K(:,nvar-z+1:nvar);
    w(k+1:k+z)=Kj'*Q(:,j);
    k=k+z;
end

% add vector w to obtain final output [vec(B);vec(Sigma);vec(w1,..,wn)]:
y=[extract_reduced_param; w];

end

