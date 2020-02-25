function [A] = makeBlockExoMatrix(info,numFor)
%input
A=zeros(info.nvar);
B=ones(numFor,info.nvar-numFor);
A((info.nvar-numFor)+1:info.nvar,1:(info.nvar-numFor))=B;
end

