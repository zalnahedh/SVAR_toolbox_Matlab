function [tt] = plotHd(HDbaseline,HDshocks,info,varSelection,shockSelection,accumP,addBaseline,colors)
%-----------------------------------------------------------------------------------------------------------
% INPUT:
% HDbaseline:      (variable,horizon,shock) /includes: lag-structure, constant and exogenous variables contribution
% HDshocks:        (variable,horizon,shock,draw) / shocks contribution
% info:            VAR information
% varSalection:    select variable for drawing historical decomposition
% shockSelection:  select shocks to plot
% accumP:          select number of periods to accumulate
%                  (To have YoY growth rates when using log(diff) transformation accumP is 4)
% addBaseline:     add baseline contribution to figure('yes' or 'no')
% colors:          set the colors for hystorical decomposition bar plot
%                  last color is for baseline and then other shocks { ...,[0 0 1]}
%------------------------------------------------------------------------------------------------------------

[v,h,s,d]=size(HDshocks);
numberOfShocks=numel(shockSelection);
shock_names=info.shocks;
number=numel(shock_names);

if strcmp(addBaseline,'yes')
  tt=zeros(h,numberOfShocks+1);  
else
  tt=zeros(h,numberOfShocks); 
end

j=0;
if strcmp(addBaseline,'yes')
    j=j+1;
    tempb=squeeze(prctile(HDbaseline(varSelection,:,:),50,3));
    tt(:,j)=tempb';
end

for i=shockSelection
    j=j+1;
    temp=squeeze(prctile(HDshocks(varSelection,:,i,:),50,4));
    tt(:,j)=temp';
end

y=info.Y(:,varSelection);

% accumulate contributions
if (accumP~=0)
    [size_1, size_2]=size(tt);
    tt_cum=zeros(size_1-accumP+1,size_2);
    y_cum=zeros(size_1-accumP+1,1);
    for i=accumP:size_1
        t1=tt(i-accumP+1:i,:);
        t2=cumsum(t1,1);
        tt_cum(i-accumP+1,:)=t2(end,:);
         y1=y(i-accumP+1:i,1);
         y2=cumsum(y1,1);
         y_cum(i-accumP+1,1)=y2(end,:);
    end
    decomp=tt_cum;
    series=y_cum;
    %date=1:size(series,1);
else
    decomp=tt;
    series=y;
    %date=1:size(series,1);   
end

%data
Xneg=decomp;
Xneg(decomp>0)=0;
Xpos=decomp;
Xpos(decomp<0)=0;

if (accumP~=0)
    time=string(info.time(accumP:end));   
else
    time=string(info.time(1:end));
end
 
x=datenum(time,'yyyyQQ');
plot(x,series,'LineWidth',3,'color','k');
datetick('x','yy','keeplimits');
hold on;

%bar plot positive side
b1=bar(x,Xpos,1,'stacked');
j=0;

if strcmp(addBaseline,'yes')
    disp(number+1)
    j=j+1;
    set(b1(j),'FaceColor',colors{number+1});
end

for i=shockSelection
    disp(i)
    j=j+1;
    set(b1(j),'FaceColor',colors{i});
end

%bar plot negative side
b2=bar(x,Xneg,1,'stacked');
j=0;

if strcmp(addBaseline,'yes')
j=j+1;
set(b2(j),'FaceColor',colors{(number+1)});
end

for i=shockSelection
    j=j+1;
    set(b2(j),'FaceColor',colors{i});
end


if strcmp(addBaseline,'yes')
   kk = cell(1,numberOfShocks+2);
else
   kk = cell(1,numberOfShocks+1);
end
kk{1}='data';
j=1;
if strcmp(addBaseline,'yes')
    j=j+1;
    kk{j}='baseline';
end
for i=shockSelection
    j=j+1;
    kk{j}=shock_names{i};
end


plot(x,series,'LineWidth',3,'color','k');
legend(kk);

hold off;
  
end

