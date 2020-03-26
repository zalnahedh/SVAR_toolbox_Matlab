function [] = saveOutputExcel(excel_name,BDraws,SigmaDraws,QDraws,WDraws,nisw,uisw,vol1,vol2,count_accept,nfinal)
%*****************************************************************************************************************
% save output of draw_nwp or draw_inwp function to excel file
%*****************************************************************************************************************
% INPUT:
%       excel_name -> excel file name
%       output of draw_nwp or draw_inwp function:
%       BDraws,SigmaDraws,QDraws, Wdraws, log_volume_1, log_volume_2, uisw and nisw
%       count_accept and nfinal
% OUTPUT:
% Excel file: 
% sheet BDraws:     BDraws
% sheet SigmaDraws: SigmaDraws
% sheet QDraws:     QDraws
% sheet Wdraws:     Wdraws
% sheet weights:    log_volume_1, log_volume_2, uisw, nisw
% counters:         count_accept, nfinal
%*******************************************************************************************************************
draws=1:size(BDraws,1);
writecell(BDraws,excel_name,'Sheet','BDraws','Range','B1');
writematrix(draws',excel_name,'Sheet','BDraws','Range','A1');

writecell(SigmaDraws,excel_name,'Sheet','SigmaDraws','Range','B1');
writematrix(draws',excel_name,'Sheet','SigmaDraws','Range','A1');

writecell(QDraws,excel_name,'Sheet','QDraws','Range','B1');
writematrix(draws',excel_name,'Sheet','QDraws','Range','A1');

writecell(WDraws,excel_name,'Sheet','WDraws','Range','B1');
writematrix(draws',excel_name,'Sheet','WDraws','Range','A1');

column_names={'log_volume_1','log_volume_2','uisw','nisw'};
writematrix(draws',excel_name,'Sheet','weigths','Range','A2');
writecell(column_names,excel_name,'Sheet','weigths','Range','B1');
writematrix(vol1,excel_name,'Sheet','weigths','Range','B2');
writematrix(vol2,excel_name,'Sheet','weigths','Range','C2');
writematrix(uisw,excel_name,'Sheet','weigths','Range','D2');
writematrix(nisw,excel_name,'Sheet','weigths','Range','E2');

counter_names={'number of accepted params','number of final params'};
counter_values={count_accept,nfinal};
writecell(counter_names',excel_name,'Sheet','counters','Range','A2');
writecell(counter_values',excel_name,'Sheet','counters','Range','B2');

end