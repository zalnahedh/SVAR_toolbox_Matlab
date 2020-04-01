function [] = saveHdExcel(excel_name,info,accumP,HDshocks,HDbaseline)
%**************************************************************************************************
% save historical decompositions to excel file (lower bound, median and upper bound)
%**************************************************************************************************
% INPUT:
%   excel_name      -> excel file name
%   info            -> VAR information
%   HDbaseline:     -> historical decomposition(variable,horizon,shock) /includes: lag-structure, constant and exogenous variables contribution
%   HDshocks:       -> historical decomposition (variable,horizon,shock,draw)/shocks contribution
% OUTPUT:
%   excel file:
%   sheet hd:      hd for each variable and for each shock 
%                   (lower bound, median, upper bound)
%**************************************************************************************************
shock_names=info.shocks;
numberOfShocks=numel(shock_names);

for i=1:size(info.names,2)
    
    varSelection=i;
    [v,h,s,d]=size(HDshocks);
    
    %**************************************************************************************************
    % median
    %*************************************************************************************************+++
    tt=zeros(h,numberOfShocks+1);  
    tempb=squeeze(prctile(HDbaseline(varSelection,:,:),50,3));
    tt(:,1)=tempb';
    j=1;
    for k=1:size(info.shocks,2)   
        j=j+1;
        temp=squeeze(prctile(HDshocks(varSelection,:,k,:),50,4));
        tt(:,j)=temp';
    end
  
    % accumulate contributions
    if (accumP(i)~=0)
        [size_1, size_2]=size(tt);
        tt_cum=zeros(size_1-accumP(i)+1,size_2);
        for k=accumP(i):size_1
            t1=tt(k-accumP(i)+1:k,:);
            t2=cumsum(t1,1);
            tt_cum(k-accumP(i)+1,:)=t2(end,:);
        end
        decomp_median=tt_cum; 
    else
        decomp_median=tt;
    end
    %**********************************************************************************************************
    % lower
    %*************************************************************************************************+++
    tt=zeros(h,numberOfShocks+1);  
    tempb=squeeze(prctile(HDbaseline(varSelection,:,:),16,3));
    tt(:,1)=tempb';
    j=1;
    for k=1:size(info.shocks,2)   
        j=j+1;
        temp=squeeze(prctile(HDshocks(varSelection,:,k,:),16,4));
        tt(:,j)=temp';
    end
  
    % accumulate contributions
    if (accumP(i)~=0)
        [size_1, size_2]=size(tt);
        tt_cum=zeros(size_1-accumP(i)+1,size_2);
        for k=accumP(i):size_1
            t1=tt(k-accumP(i)+1:k,:);
            t2=cumsum(t1,1);
            tt_cum(k-accumP(i)+1,:)=t2(end,:);
        end
        decomp_lower=tt_cum; 
    else
        decomp_lower=tt;
    end
    %*********************************************************************************************************
    %upper
    %*************************************************************************************************+++
    tt=zeros(h,numberOfShocks+1);  
    tempb=squeeze(prctile(HDbaseline(varSelection,:,:),84,3));
    tt(:,1)=tempb';
    j=1;
    for k=1:size(info.shocks,2)   
        j=j+1;
        temp=squeeze(prctile(HDshocks(varSelection,:,k,:),84,4));
        tt(:,j)=temp';
    end
  
    % accumulate contributions
    if (accumP(i)~=0)
        [size_1, size_2]=size(tt);
        tt_cum=zeros(size_1-accumP(i)+1,size_2);
        for k=accumP(i):size_1
            t1=tt(k-accumP(i)+1:k,:);
            t2=cumsum(t1,1);
            tt_cum(k-accumP(i)+1,:)=t2(end,:);
        end
        decomp_upper=tt_cum; 
    else
        decomp_upper=tt;
    end
 %*********************************************************************************************************
 % time
 %*********************************************************************************************************
    if (accumP(i)~=0)
        time=string(info.time(accumP(i):end));   
    else
        time=string(info.time(1:end));
    end
 %*********************************************************************************************************
    
    for j=1:size(decomp_median,2) 
        
        decomp=[decomp_lower(:,j) decomp_median(:,j) decomp_upper(:,j)];
        
        var_name=strcat('variable: ',info.names{i});
        if (j==1)
           shock_name=strcat('contribution/baseline'); 
        else
            shock_name=strcat('contribution/shock/ ',info.shocks{j-1});
        end
        [shock_position,var_position]=excelPositionIrf(info,i,j);
        shock_position_letter=xlcolumnletter(shock_position);
        info_position=strcat(shock_position_letter,num2str(var_position));
        column_position=strcat(shock_position_letter,num2str(var_position+1));
        matrix_position=strcat(shock_position_letter,num2str(var_position+2));
        column_names={'low','median','up'};
        info_names={var_name,shock_name};   
        writecell(info_names,excel_name,'Sheet','hd','Range',info_position);
        writecell(column_names,excel_name,'Sheet','hd','Range',column_position);
        writematrix(decomp,excel_name,'Sheet','hd','Range',matrix_position);
    end
        
        
        
    
 
end



end %end
