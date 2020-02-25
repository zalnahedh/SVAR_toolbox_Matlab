function [ Y,X ] = dataForVar(info)
    if nargin~=1
        error(strcat('Wrong number of arguments for function DataForVar. You put: ',num2str(nargin),'and 1 argument is required'));
    end
    p=info.nlag;
    nex=info.nex;
    endo=info.endo;
    if nex~=0
        exo=info.exo;
    end
    cte=info.cte;
  
  
    % put lagged data into matrix X 
    X=lag0(endo,1); 
    for i=2:p
        X=[X lag0(endo,i)];
    end
    
    % if there are exogenous variables put them into matrix X
    if nex~=0
        X=[X exo];
    end
    
    % if there is constant put it into matrix X
    if cte~=0
        Xones=ones(size(endo,1),1);
        X=[X Xones];
    end
    
    % ignore first p rows of zeros and create matrices X and Y
    X=X(p+1:end,:);
    Y=endo(p+1:end,:);
   
end

