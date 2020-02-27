function output = F_map(x,info)
%------------------------------------------------------------------------
% input: x=[vec(A0);vec(Aplus)] structural representation
% output: matrix F
% F is evaluated at point x
%----------------------------------------------------------------------

nvar=info.nvar;
restrictions=info.restrictions;

n_squared=nvar*nvar;

m=info.m;

%initialize matrix F
F=[];

% Add A0 to F if cell restrictions contains 'A0'
if any(strcmp(restrictions,'A0'))
        % extract matrix A0 from vectorized implementation x
        A0=reshape(x(1:n_squared),nvar,nvar);
        % add matrix AO to matrix F
        F=[F;A0];
end    

% Add Aplus to F if cell restrictions contains 'Aplus'
if any(strcmp(restrictions,'Aplus'))
        % extract matrix Aplus from vectorized implementation x
        Aplus=reshape(x(n_squared+1:end),m,nvar);
        F=[F;Aplus];
end     


% calculate impulse responses for given structural representation x
irfs=irf_short(x,info);

% calculate F
num_res=size(restrictions,2);
for i=1:num_res
    s=restrictions{i};
    if ((~strcmp(s,'A0')) && (~strcmp(s,'Aplus')))
        F=[F;irfs(s)];
    end
end

output =F;

end
