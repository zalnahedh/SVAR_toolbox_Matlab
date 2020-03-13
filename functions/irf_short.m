function output = irf_short(x,info)


nvar=info.nvar;
p=info.nlag;
restrictions=info.restrictions;
num_res=size(restrictions,2);
 
n_squared=nvar*nvar;
nx=numel(x); 
m=nx/nvar-nvar;

A0=reshape(x(1:n_squared),nvar,nvar);

Aplus=reshape(x(n_squared+1:end),m,nvar);
 
A=cell(p,1);
for i=1:p
   
    A{i}=Aplus((i-1)*nvar+1:i*nvar,:);
end
 

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

for i=1:num_res
    s=restrictions{i};
    if ((~strcmp(s,'A0')) && (~strcmp(s,'Aplus')))
        if (strcmp(s,'0'))
           L(s)=A0_inv'; 
        end
        if ((strcmp(s,'long_run')))
           L(s)=A0';
           for j=1:p
               L(s)=L(s)-A{j}';
           end
           %L(s)=inv(L(s)); 
           L(s)=L(s)\eye(info.nvar);
        end
        if ((~strcmp(s,'0')) && (~strcmp(s,'long_run')))
            ss=str2num(s);
            L(s)=(A0_inv*J'*F^(ss)*J)';
        end      
    end
end
 


output=L;
 
end
