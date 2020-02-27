function output=ZF(x,info)
%------------------------------------------------------------------------
% input: x=[vec(A0);vec(Aplus)] structural representation
% output: ZF=cell(n,1)
% ZF is product of Z matrix nad F evaluated at given structural rep. x
%----------------------------------------------------------------------

Z=info.Z;
nvar=info.nvar;

% calculate matrix F at point x
F=F_map(x,info);

% calculate product of Z and F called ZF
ZF=cell(nvar,1);
for i=1:nvar
   ZF{i}=Z{i}*F;
end

output=ZF;

end
