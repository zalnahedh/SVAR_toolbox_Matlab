function s=stand(Y,X,info)
    if nargin~=3
        error('Pogresan broj argumenata');
    end
    p=info.nlag; %number of lags
    n=info.nvar; %number of variables
   
    s=[];
    for i=1:n
        y=Y(:,i);
        x=X(:,i); 
        Xones=ones(size(x,1),1);
        x=[x Xones];
        b0=inv(x'*x)*(x'*y);
        s_pom=sqrt(((y-x*b0)'*(y-x*b0))/(size(y,1)-2));
        s=[s s_pom];
    end
end