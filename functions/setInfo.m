function info = setInfo(nvar,m,nlag,nex,Z,S,cte,endo_variables,exo_variables,variable_names,shock_names,impact);
 
%overall number of zero restrictions
numZeroRes=0;
%overall number of sign restrictions
numSignRes=0;
% W_j is some fixed (nvar+1-j-z_j) x n matrix
W_matrix=cell(nvar,1);

% draw w_i from the uniform distribution over the unit sphere in
% R^{nvar+1-j-z_j} for each 1 <= j <= nvar
% compute the dimension of space in which the spheres live that map into
% the orthogonal matrices satisfying the restrictions
dimSphere=0;

for j=1:nvar
    %number of rows of matrix Z{j} that represents zero restrictions for
    %j-th structural shock
    z_j=size(Z{j},1); 
    W_matrix{j}=randn(nvar+1-j-z_j,nvar);
    numZeroRes=numZeroRes+z_j;
    dimSphere=dimSphere+(nvar+1-j-z_j);
end

for j=1:nvar
   s_j=size(S{j},1);
   numSignRes=numSignRes+s_j;
end
  
% number of endogenous variables
info.nvar=nvar; 
% number of exogenous variables
info.nex =nex;
% lag 
info.nlag=nlag;
% constant inclued(1) or not included (0)
info.cte =cte;
% dimension of matrix Aplus in structural representation is (m x nvar)
info.m=m;  
% zero restrictions selection matrix
info.Z=Z;
% sign restrictions selection matrix
info.S=S;
% dimension of matrix C in SVAR(p) notation is (k x nvar), last column is column of ones if constant is included
info.k=nex+cte; 

info.dimSphere=dimSphere;
info.numZeroRes=numZeroRes;
info.numSignRes=numSignRes;
info.W_matrix=W_matrix;

% handle exogenous variables if exist
if nex~=0  
    % matrix of exogenous variables (full length after transformation(diff))
    info.exo        =exo_variables; 
end

% matrix of endogenous variables (full length after transformation(diff))
info.endo=endo_variables;
% cell array of imposed restrictions (periods: 0,1,..;'long_run'; 'A0' or 'Aplus')
info.restrictions=keys(impact);
% cell array of endogenous variables names
info.names=variable_names;
% structural shocks names (defined by imposing sing and zero restrictions)
info.shocks=shock_names;
end

