function [hd_a,hd_b,hd_c] = hd(x,irf_matrix,info )


nvar=info.nvar;
p=info.nlag;
horizons=info.horizons;
nex=info.nex;
cte=info.cte;
 
n_squared=nvar*nvar;
nx=numel(x); 
m=nx/nvar-nvar;
k=nex+cte;
endo=info.endo;



% extract matrix A0 from vectorized implementation x
A0=reshape(x(1:n_squared),nvar,nvar);
% extract matrix Aplus from vectorized implementation x
Aplus=reshape(x(n_squared+1:end),m,nvar);
 
% extract matrices A(1),...,A(p) from matrix Aplus
A=cell(p,1);
for i=1:p
    % extract submatrix A{i} from matrix Aplus
    A{i}=Aplus((i-1)*nvar+1:i*nvar,:);
end

%extract matrix C from matrix Aplus
C=Aplus(p*nvar+1:end,:);

[Y,X]   =dataForVar(info);
errors=Y*A0-X*Aplus;
zeta=X(:,nvar*p+1:end);

%calculate F matrix and J matrix 
identityMatrix1=eye(nvar,nvar);
identityMatrix2=eye(nvar*(p-1),nvar*(p-1));
zeroMatrix=zeros(nvar*(p-1),nvar);
if p==1
    J=identityMatrix1;
    A0_inv=inv(A0);
    F=A{1}*A0_inv;
else
    J=[identityMatrix1;zeroMatrix];
    F=zeros(nvar*p,nvar*p);
    A0_inv=inv(A0);
    for i=1:p
       F(nvar*(i-1)+1:i*nvar,1:nvar)=A{i}*A0_inv;
    end
    F(1:(p-1)*nvar,nvar+1:end)=identityMatrix2;
end

hd_a=zeros(nvar,horizons);
hd_b=zeros(nvar,horizons,k);
hd_c=zeros(nvar,horizons,nvar);

y0=endo(p,:)';
for i=p-1:-1:1
    y0=[y0;endo(i,:)'];
end

help=zeros(nvar,1);
for t=0:(horizons-1)
        help=J'*(F')^(t+1)*y0;
        for kk=1:nvar
            hd_a(kk,t+1)=help(kk,1);
        end
end

for kk=1:nvar
    for t=1:(horizons)
        for j=1:nvar
            suma=0;
            for i=0:(t-1)
                r=t-i;
                new=irf_matrix(kk,j,i+1)*errors(r,j);
                suma=suma+new;
            end
            hd_c(kk,t,j)=suma;
        end
    end
end

help=zeros(nvar,k,horizons);
for t=0:(horizons-1)
    help(:,:,t+1)=irf_matrix(:,:,t+1)*C';
end

for kk=1:nvar
    for t=1:(horizons)
        for j=1:k
            suma=0;
            for h=0:(t-1)
                r=t-h;
                new=help(kk,j,h+1)*zeta(r,j);
                suma=suma+new;
            end
            hd_b(kk,t,j)=suma;
        end
    end
end


end

