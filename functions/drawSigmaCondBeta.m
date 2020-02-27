function [sigma] = drawSigmaCondBeta(info,beta,scaleMatrix0,alfa0)

k=info.nvar*info.nlag+info.nex+info.cte; %total number of VAR coefficients per each equation
e=info.Y-info.X*reshape(beta,k,info.nvar);

scaleMatrix=e'*e+scaleMatrix0;
sigma=iwpq(info.T+alfa0,inv(scaleMatrix));
end

