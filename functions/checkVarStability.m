function [stability] = checkVarStability(info,beta)
stability=1; %var is stable
n=info.nvar;
p=info.nlag;
FF=zeros(n*p,n*p);
FF(n+1:n*p,1:n*(p-1))=eye(n*(p-1),n*(p-1));
temp=reshape(beta,n*p+info.nex+info.cte,n);
temp=temp(1:n*p,1:n)';
FF(1:n,1:n*p)=temp;
ee=max(abs(eig(FF)));
if ee>1
    % var is not stable
    stability=0;
end

end

