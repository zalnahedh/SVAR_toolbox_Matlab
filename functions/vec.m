function output=vec(Y)

% --Description---------------------------------------------------------------
% Function vec transforms matrix Y into vector output taking column by
% column of matrix Y 
%--------------------------------------------------------------------------
output=[];
for i=1:size(Y,2)
    output=[output;Y(:,i)];
end
end

