function output=zeroRestrictions(x, info)
%
%
%
nvar=info.nvar;
numZeroRes=info.numZeroRes;
ZF=info.ZF(x,info);

output=zeros(numZeroRes,1);
k=1;
for i=1:nvar
    z=size(ZF{i},1);
    output(k:k+z-1)=ZF{i}(:,i);
    k=k+z;
end
end