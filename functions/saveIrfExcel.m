function [] = saveIrfExcel(excel_name,info,IrfDraws)
%**************************************************************************************************
% save irfs to excel file (lower bound, median and upper bound)
%**************************************************************************************************
% INPUT:
%   excel_name      -> excel file name
%   info            -> VAR information
%   IrfDraws        -> impulse response function
% OUTPUT:
%   excel file:
%   sheet irf:      irf for each variable and for each shock 
%                   (lower bound, median, upper bound)
%**************************************************************************************************+
for j=1:size(info.shocks,2)   
    for i=1:size(info.names,2)
        var_name=strcat('variable: ',info.names{i});
        shock_name=strcat('shock: ',info.shocks{j});
        [shock_position,var_position]=excelPositionIrf(info,i,j);
        shock_position_letter=xlcolumnletter(shock_position);
        info_position=strcat(shock_position_letter,num2str(var_position));
        column_position=strcat(shock_position_letter,num2str(var_position+1));
        matrix_position=strcat(shock_position_letter,num2str(var_position+2));
        A=IrfDraws;
        A(i,:,:,:)=cumsum(IrfDraws(i,:,:,:),3);
        D=reshape(A(i,j,1:end,:),info.horizons,info.nfinal);
        temp=[prctile(D',[16 50 84])];
        column_names={'low','median','up'};
        info_names={var_name,shock_name};   
        writecell(info_names,excel_name,'Sheet','irf','Range',info_position);
        writecell(column_names,excel_name,'Sheet','irf','Range',column_position);
        writematrix(temp',excel_name,'Sheet','irf','Range',matrix_position);
    end
end

end