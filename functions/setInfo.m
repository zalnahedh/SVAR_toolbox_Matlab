function info = setInfo(info)
%overall number of zero restrictions
numZeroRes=0;
%overall number of sign restrictions
numSignRes=0;
% W_j is some fixed (nvar+1-j-z_j) x n matrix
W_matrix=cell(info.nvar,1);

% draw w_i from the uniform distribution over the unit sphere in
% R^{nvar+1-j-z_j} for each 1 <= j <= nvar
% compute the dimension of space in which the spheres live that map into
% the orthogonal matrices satisfying the restrictions
dimSphere=0;

for j=1:info.nvar
    %number of rows of matrix Z{j} that represents zero restrictions for
    %j-th structural shock
    z_j=size(info.Z{j},1); 
    W_matrix{j}=randn(info.nvar+1-j-z_j,info.nvar);
    numZeroRes=numZeroRes+z_j;
    dimSphere=dimSphere+(info.nvar+1-j-z_j);
end

for j=1:info.nvar
   s_j=size(info.S{j},1);
   numSignRes=numSignRes+s_j;
end

info.dimSphere=dimSphere;
info.numZeroRes=numZeroRes;
info.numSignRes=numSignRes;
info.W_matrix=W_matrix;

end

