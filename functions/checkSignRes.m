function y =checkSignRes(x,info)
%------------------------------------------------------------------------
% input: x=[vec(A0);vec(Aplus)] structural representation
% output:  True if sign restrictions satisfied
%          False if sign restrictions are not satisfied
%----------------------------------------------------------------------

S=info.S;
numSignRes=info.numSignRes;
nvar=info.nvar;

% calculate matrix F at point x
F=F_map(x,info);

help=zeros(numSignRes,1);

ib=1;
for j=1:nvar
    ie=ib+size(S{j},1);  
    help(ib:ie-1)=S{j}*F(:,j);
    ib=ie;
end

if (sum(help>0))==numSignRes
    y=true;
else
    y=false;
end
    