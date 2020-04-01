function [] = saveInfoExcel(excel_name,info)
%*******************************************************************************************************
% saving to excel file 
%******************************************************************************************************
% INPUT:
%   excel_name      -> excel file name
%   info            -> VAR information
%OUTPUT:
%  Sheet 'info' :
%   variable's names, shock's names, variable's names, restrictions
%   'maxDraws','maxQ','finalDraws','iter_show','index','nlag','cte','nvar',
%   'nex','k','m','dimSphere','numZeroRes','numSignRes','horizons'
%  Sheet 'data' :
%    Transformed and truncated endogenous variables
%  Sheet 'restrictions' :
%    sign and zero restrictions matrices
%*********************************************************************************************************
variable_names={' variable names: '};
writecell(variable_names,excel_name,'Sheet','info','Range','A2');
writecell(info.names,excel_name,'Sheet','info','Range','B2');

% save shock's names
shock_names={'shock names: '};
writecell(shock_names,excel_name,'Sheet','info','Range','A3');
writecell(info.shocks,excel_name,'Sheet','info','Range','B3');

% save restrictions
restriction_names={'restrictions: '};
writecell(restriction_names,excel_name,'Sheet','info','Range','A4');
writecell(info.restrictions,excel_name,'Sheet','info','Range','B4');

% save other settings
info_names={'maxDraws','maxQ','finalDraws','iter_show','index','nlag','cte','nvar','nex','k','m','dimSphere','numZeroRes','numSignRes','horizons'};
info_values={info.maxDraws,info.maxQ,info.finalDraws,info.iter_show,info.index,info.nlag,info.cte,info.nvar,info.nex,info.k,info.m,info.dimSphere,info.numZeroRes,info.numSignRes,info.horizons};
writecell(info_names',excel_name,'Sheet','info','Range','A5');
writecell(info_values',excel_name,'Sheet','info','Range','B5');

% save data for VAR (transformed and truncated)
writecell({'transformed and truncated (by number of lags'},excel_name,'Sheet','endogenous data','Range','B1');
writecell(info.names,excel_name,'Sheet','endogenous data','Range','B2');
writematrix(info.Y,excel_name,'Sheet','endogenous data','Range','B3');
writecell(info.time,excel_name,'Sheet','endogenous data','Range','A3');

% save sign and zero restrictions impact matrices
restrictions=keys(info.impact);
k=1;
for i=1:size(info.impact,1)
    s=restrictions{i};
    output_matrix=info.impact(s);
    text_position=strcat(xlcolumnletter(k),'1');
    matrix_position=strcat(xlcolumnletter(k+1),'1');
    text = strcat('restriction: ',s);
    writecell({text},excel_name,'Sheet','restrictions','Range',text_position);
    writematrix(output_matrix,excel_name,'Sheet','restrictions','Range',matrix_position);
    k=k+size(output_matrix,2)+1;
end
end
