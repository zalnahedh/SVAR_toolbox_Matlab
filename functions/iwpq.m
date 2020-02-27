function output=iwpq(d,S)
%------Decription----------------------------------------------------
% simulate matrix from IW distribution with d degrees of freedom and scale
% matrix S
%----------------------------------------------------------------
n=size(S,1);
z=zeros(d,n);
for i=1:d
    help=randn(n,1);
    z(i,:)=(chol(S)'*help)';
end
output=inv(z'*z);
end

