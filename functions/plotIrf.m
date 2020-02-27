function [] = plotIrf(irfArray,info,shockList,variableList,cumulativeList)

% irfArray:   calculated IRFs (variable, shock, horizont, draw)
% info:         VAR info
% shockList:    list of shocks, for example: 1:4 or [1,3]
% variableList: list of variables, for example: 1:6, [1,4,5]
% cumulativeList: list  of variables for which the corresponding impulse
% responses need to be accumulated

A=irfArray;

for k=cumulativeList
    A(k,:,:,:)=cumsum(irfArray(k,:,:,:),3);
end

br=0;
nfinal=info.nfinal;
horizon=info.horizons;
names=info.names;
shocks=info.shocks;

numberOfShocks=numel(shockList);
numberOfVariables=numel(variableList);

for j=shockList
    for i=variableList
        br=br+1;
        subplot(numberOfShocks,numberOfVariables,br)
        D=reshape(A(i,j,1:end,:),horizon,nfinal);
        temp=[prctile(D',[50 16 84])];
        plot(temp(1,:),'LineWidth',1.5)
        jbfill(1:horizon,temp(2,:),temp(3,:),[0.1 0.6 0.9],[1 0 0],1,0.3);
        title([shocks{j},' ==> ', names{i}]);
        grid on
    end
end

end

