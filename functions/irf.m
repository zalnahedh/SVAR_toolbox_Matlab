function output = irf(x,info)


nvar=info.nvar;
p=info.nlag;
horizons=info.horizons;
 
n_squared=nvar*nvar;
nx=numel(x); 
m=nx/nvar-nvar;

A0=reshape(x(1:n_squared),nvar,nvar);
Aplus=reshape(x(n_squared+1:end),m,nvar);
 

A=cell(p,1);
for i=1:p
    A{i}=Aplus((i-1)*nvar+1:i*nvar,:);
end
 
%extract matrix C from matrix Aplus
%C=Aplus(p*nvar+1:end,:);
 
%IRF representation
L=containers.Map();
 
identityMatrix1=eye(nvar,nvar);
identityMatrix2=eye(nvar*(p-1),nvar*(p-1));
zeroMatrix=zeros(nvar*(p-1),nvar);

if p==1
    J=identityMatrix1;
    %A0_inv=inv(A0);
    A0_inv=A0\eye(info.nvar);
    F=A{1}*A0_inv;
else
    J=[identityMatrix1;zeroMatrix];
    F=zeros(nvar*p,nvar*p);
    %A0_inv=inv(A0);
    A0_inv=A0\eye(info.nvar);
    for i=1:p
       F(nvar*(i-1)+1:i*nvar,1:nvar)=A{i}*A0_inv;
    end
    F(1:(p-1)*nvar,nvar+1:end)=identityMatrix2;
end

% calculate impulse responses:
for i=1:horizons 
    L(int2str(i))=(A0_inv*J'*F^i*J)';
end

L('0')=A0_inv';

L('long_run')=A0';
for i=1:p
    L('long_run')=L('long_run')-A{i}';
end
%L('long_run')=inv(L('long_run'));
 L('long_run')=(L('long_run'))\eye(info.nvar);
output=L;
 
end
