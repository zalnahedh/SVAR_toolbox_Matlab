
function output = h_tilda(X, info)
%
% X - n x n matrix
%
% output- n x n matrix such that Y'*Y = 0.5*(X+X')
%
 
output=info.h(0.5*(X+X'));
end
