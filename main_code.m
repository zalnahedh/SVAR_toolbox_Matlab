
%% BAYESIAN VECTOR AUTOREGRESSION MODELS (structural analysis )
%  BVAR with Block Exogeneity,Sign and Zero Restrictions
% *****************************************************************************************************************
%
% Author : ANTE COBANOV (ante.cobanov@yahoo.com)
% GitHub:  acobanov
%
%*******************************************************************************************************************
% Comments: 
% Importance sampling algorithm for imposing sign and zero restrictions is implemented as in Jonas E.Arias,
% Juan F.Rubio-Ramirez and Daniel F.Waggoner(2018): "Inference Based on SVARs Identified with Sign and Zero 
% Restrictions: Theory and Applications",Econometrica, March 2018, Volume
% 86. I am very grateful to my colleagues: Karlo Kotarac & Davor Kunovac
% for the great suggestions,ideas, comments and helping me to fix some programming
% bugs.

% UNDER CONSTRUCTION:
% Problems with large values of importance sampling weights!!
% Add tests: "Testing the assumptions behind importance sampling" 
%             Koopman, Siem Jan; Shepard; Neil; Creal; Drew (Journal of
%             Econometrics)

% Modifications needed: "Inferance in Bayesian Proxy-SVARs"
%                       Arias; Rubio-Ramirez; Waggoner (2019)
%***************************************************************************************************************************
%
% Matlab code includes importance sampling algorithm under Block Exogeneity,Sign and Zero Restrictions with possibility to include 
% exogenous variables and to extend restrictions on any admissible function of structural parameters
%
%******************************************************************************************************************************
%% DATA AND WORKING DIRECTORY
%=========================================================================================================================================================
clear variables;
close all;
userpath('clear');
clc;
%tic;

rng('default'); % reinitialize the random number generator to its startup configuration
rng(0);         % set seed

currdir=pwd;
%cd ..
get_help_dir_currdir=pwd;
addpath([get_help_dir_currdir,'/functions']); % set path to helper functions
cd(currdir)

% reading data from excel
input_file='data/data.xlsx';
info.format='yyyyQQ';

sheets={'hr'};

for sheet_index=1:numel(sheets)
    
    sheet=sheets{sheet_index};
    [data,txt]=xlsread(input_file,sheet);
    names=txt(1,2:end);
    time=txt(2:end,1);
    data=data(:,1:end);
    date=datenum(time,info.format); 

    % selecting endogenous variables and their order
    % to select : second, first, fourth and third variable in excel put [2,1,4,3]
    % order is important for imposing block-exogeneity (put first "small economy" variables)
    endo_variables =data(1:end,1:4);

    % cell array of endogenous variables names
    info.names=names(1,1:4);
    % import exogenous variables (put empty matrix [] if there are no exogenous variables)
    exo_variables= []; 

    % transform endogenous variables
    % matrix of endogenous variables (full length after transformation(diff))
    endo_variables=diff(log(endo_variables));
    %transform exogenous variables if needed!!!
    %exo_variables=diff(log(exo_variables));

    start =size(data,1)-size(endo_variables,1);
    date=date(start+1:end,:);
    time=time(start+1:end,:);

    %========================================================================================================================================================
    %% MODEL SETUP: 
    %========================================================================================================================================================
    %  Notation: Write data in Rubio, Waggoner, and Zha (RES 2010)'s notation with modification: vector of exogenous variables z(t)  
    %  Label:  't' in xt,yt,et,zt,Aplust,etc means transponse
    %  Label: '(t)' in xt(t),yt(t),et(t),zt(t),etc means time
    %
    % SVAR(p):                    yt(t) A0 = yt(t-1) A1 + ... +yt(t-nlag) Anlag +zt(t) C + et(t)  for t=1...,T;
    % structural rep.(AO,Aplus):  yt(t) A0 = xt(t) Aplus + et(t)  for t=1...,T;
    % reduced-form rep.(B,Sigma): yt(t)= xt(t) B + ut(t)          for t=1...,T;
    % where:                      xt(t) = [yt(t-1), ... , yt(t-nlag),zt(t)];
    %                             Aplust(t)=[A1t, ... , Anlagt, Ct]
    %                             B=Aplus(A0)^(-1) &  Sigma=(A0 A0t)^(-1)

    % matrix notation:            Y = X*B + U;
    % where:                      Y= [yt(nlag+1);yt(nlag+2); ... ;yt(end)]
    %----------------------------------------------------------------------------------------------------------------------------------------------------
    prior                          =  "inwp";    % set inwp (independent NW) or nwp (NW)
    info.maxDraws                  =  1000000;  % maximum number of importance sampling orthogonal-reduced-form (B,Sigma,Q) draws
    info.maxQ                      =  100;      % draw several matrices Q for the same (B,Sigma) to speed up proccess
    info.finalDraws                =  200;      % effective size of importance sampling draws (smaller than accepted draws)
    info.iter_show                 =  50;       % display iteration every iter_show draws
    info.index                     =  40;       % define  horizons for IRFs 
    info.nlag                      =  2;        % number of lags in a SVAR model
    info.cte                       =  1;        % set equal to 1 if a constant is included; 0 otherwise
    %-----------------------------------------------------------------------------------------
    info.nvar      = size(endo_variables,2);           % number of endogenous variables
    info.nex       = size(exo_variables,2);            % number of exogenous variables in a SVAR model
    info.k         = info.nex + info.cte;              % dimension of matrix C in SVAR(p) notation is (k x nvar), last column is column of ones if constant is included
    info.m         = info.nvar*info.nlag + info.k;     % dimension of matrix Aplus in structural representation is (m x nvar)
    %----------------------------------------------------------------------------------------------------------------------------------------------------

    %=====================================================================================================================================================
    %% DEFINE ZERO AND SIGN RESTRICTIONS:
    %=====================================================================================================================================================
    %  Define impact patern matrix (impact) and selection matrices: zero (Z) and sign (S) 
    %  To define impact matrix you can define patern for: 'irf', 'structural' or 'both' 
    %  If 'irf', user should define restrictions to any horizon of IRFs, use strings : '0', '1', '2', etc or 'long_run' 
    %  If 'structural', user should define restrictions to structural matrices,use strings : 'A0','Aplus'
    %------------------------------------------------------------------------------------------------------------------------------------------------

    info.type_of_restrictions = 'structural';    %set the type of restrictions: 'irf', 'structural' or 'both'

    % structural shocks names (defined by imposing sing and zero restrictions)
    info.shocks={'HR demand','HR supply','EA demand', 'EA supply'}; 

    %set the restrictions ({'0','long_run','A0','Aplus'})
    % initialize impact pattern
    % define restrictions (in matrix form) for each element in vector of restrictions
    impact=containers.Map();
    impact('0')=[1 1 0 0       % 1: HR aggregate demand   
                 1 -1 0 0      % 2: HR aggregate supply
                 9 9 1 1       % 3: EA aggregate demand
                 9 9 1 -1 ];   % 4: EA aggregate supply 
    impact('long_run')=[0 9 9 9       
                         9 9 9 9      
                         9 9 0 9      
                         9 9 9 9 ] ;   
    % calculate selection matrices: ZZ (zero) and SS (sign)
    [SS,ZZ]=get_SZ(impact);
    % cell array of imposed restrictions (periods: 0,1,..;'long_run'; 'A0' or 'Aplus')
    % Take a look at the F_map function, restrictions in vector of restrictions need  to be in this order: A0, Aplus, IRF
    info.restrictions=keys(impact);
    % zero restrictions selection matrix
    info.Z=ZZ;
    % sign restrictions selection matrix
    info.S=SS;

    % handle exogenous variables if exist
    if info.nex~=0  
        % matrix of exogenous variables (full length after transformation(diff))
        info.exo        =exo_variables; 
    end
    info.endo =endo_variables;
    %=========================================================================================================================================
    %% SETTING UP INFO:
    %=========================================================================================================================================
    % this helps to organize program code better 
    info=setInfo(info);

    % prepare the DATA
    [Y,X]   =dataForVar(info);
    T       =size(Y,1);
    info.T=T;
    info.time=time(info.nlag+1:end);
    info.date=date(info.nlag+1:end);
    info.Y=Y;
    info.X=X;

    % define  horizons for HD 
    info.horizons=T;

    % help functions (some to compute importance sampler weights(isw))
    info.h                  =@(x)chol(x);
    info.ZF                 =@(x,y)ZF(x,y);
    info.strToW             =@(x)gf_h_map(x,info);
    info.zeroRestr          =@(x)zeroRestrictions(x,info);
    info.strToRed           =@(x)f_h_map(x,info);
    info.strToIrf           =@(x)structToIrfRecursive(x,info);
    info.irfToStr           =@(x)irfToStructRecursive(x,info);
    %===========================================================================================================================================
    %% MINNESOTA PRIOR (use Independent Normal Wishart prior)
    %********************************************************************************************************************************************
    if strcmp(prior,'inwp')       
            lambda.lambda1  =0.1;          % hyperparameter for covariance matrix of reduced-form coefficients
            lambda.lambda2  =0.5;          % hyperparameter for covariance matrix of reduced-form coefficients
            lambda.lambda3  =1;            % hyperparameter for covariance matrix of reduced-form coefficients
            lambda.lambda4  =100;          % hyperparameter for covariance matrix of reduced-form coefficients
            lambda.lambda5  =0.001;        % hyperparameter for block-exogeneity
            lambda.firstLagCoef=0.8;       % mean prior ; put 0.8 for stationary variables 
            lambda.stand =stand(Y,X,info); % calculate OLS residual variance of the auto-regressive models estimated for each variable
     end
    %***********************************************************************************************************************************
    %% BLOCK EXOGENEITY (use Independent Normal Wishart prior)
    %=======================================================================================================================================
    if strcmp(prior,'inwp') 

            blockExo.bex=1;         % put 1 if there are block-exogeneity restrictions, otherwise put 0 
            blockExo.option=2;      % put 1 if you want to impose block-exogeneity by number of foreign variables
                                    % put 2 if you want to impose  by user defined restrictions
                                    % option 2 allows to have more blocks
            numForeignVar=2;        % number of foreign variables

            % initialize block-exogeneity (no restrictions imposed)
            blockExo.bexMatrix=zeros(info.nvar);

            % option 1 (imposing block-exogeneity by number of foreign varables)
            if ((blockExo.option==1)&&(blockExo.bex==1))
                blockExo.bexMatrix=makeBlockExoMatrix(info,numForeignVar);
            end

            % option 2 (imposing user-defined )
            if ((blockExo.option==2)&&(blockExo.bex==1))
                % create your block-exogeneity matrix here 
                % put 1 to impose block-exogeneity restrictions(1 means "zero" i.e. no response)
                blockExo.bexMatrix = [0 0 1 1         
                                      0 0 1 1     
                                      0 0 0 0     
                                      0 0 0 0];
                blockExo.bexMatrix=blockExo.bexMatrix';
            end  

            %decide number of draws to be rejected until convergence
            d2.burn=1000;

            % impose prior (Minnesota and block-exogeneity)
            [d2.beta0,d2.omega0] = minnesotaPrior(info,lambda,blockExo);

            % alternatively put zero mean prior 
            %qq=info.nvar*(info.nvar*info.nlag+info.nex+info.cte);
            %beta0=zeros(q,1);

            % alternatively put OLS estimate for mean prior
            %beta0=vec(invpd((info.X)'*info.X)*(info.X)'*(info.Y));

            d2.alfa0=info.nvar+2;
            % choose scaleMatrix0 
            %scaleMatrix0=eye(info.nvar);
            d2.scaleMatrix0=(d2.alfa0-info.nvar-1)*diag(lambda.stand.^2);   
    end
    %=====================================================================================================================================================
    %% SET PRIOR AND POSTERIOR  ( NORMAL WISHART PRIOR )
    %=====================================================================================================================================================
    if strcmp(prior,'nwp')   
            %  Prior distribution over the reduced-form parameters is normal-inverse-Wishart (NIW)
            %  "Bar" noration is used for priors and Inv for inverse!
            %  NIW is characterized by 4 parameters: phi matrix,scalar nu, psi matrix, omega matrix 

            % phi matrix:   symmetric and positive definite  (nvar x nvar) (scale matrix of IW-distribution)
            % scalar nu:    degree of freedom of IW-distribution (>=n)
            % psi matrix:       (m x nvar) (conditionally normal)
            % omega matrix:     symmetric and positive definite (m x m) (conditionally normal)

            d1.phiBar       = zeros(info.nvar);
            d1.nuBar        = 0;     
            d1.psiBar       = zeros(info.m,info.nvar);
            d1.omegaBarInv  = zeros(info.m);

            % If the prior is NIW then the posterior distribution over the reduced-form
            % parameteres is also NIW with the given transformation formulas (notation "Tilda"):
            d1.nuTilda         = info.T +d1.nuBar;
            d1.omegaTilda      = (X'*X  + d1.omegaBarInv)\eye(info.m);
            d1.omegaTildaInv   =  X'*X  + d1.omegaBarInv;
            d1.psiTilda        = d1.omegaTilda*(X'*Y + d1.omegaBarInv*d1.psiBar);
            d1.phiTilda        = Y'*Y + d1.phiBar + d1.psiBar'*d1.omegaBarInv*d1.psiBar - d1.psiTilda'*d1.omegaTildaInv*d1.psiTilda;
            d1.phiTilda        = (d1.phiTilda'+d1.phiTilda)*0.5; %symmetric

            % definition to facilitate the draws from B|Sigma
            d1.cholOmegaTilda  = (info.h(d1.omegaTilda))'; % this matrix is used to draw B|Sigma below
    end
    %=====================================================================================================================================================
    %% NORMAL WISHART PRIOR  (algorithm)
    %**************************************************************************************************************************
    if strcmp(prior,'nwp')
      [BDraws,SigmaDraws,QDraws,WDraws,nisw,uisw,vol1,vol2,count_accept,nfinal] = draw_nwp(info,d1);
    end

    %**************************************************************************************************************************
    %% INDEPENDENT NORMAL WISHART PRIOR  (algorithm)
    %****************************************************************************************************************
    if strcmp(prior,'inwp')
      [BDraws,SigmaDraws,QDraws,WDraws,nisw,uisw,vol1,vol2,count_accept,nfinal] = draw_inwp(info,d2);
    end

    %**********************************************************************************************************
    %% COMPUTE impulse response functions and hystorical decompositions
    %****************************************************************************************************************

    info.nfinal=nfinal;

    %initialization
    IrfDraws=zeros(info.nvar,info.nvar,info.horizons,info.nfinal); % impulse response function
    HDaDraws=zeros(info.nvar,info.horizons,info.nfinal);           % lag-structure (hd)
    HDbDraws=zeros(info.nvar,info.horizons,info.k,info.nfinal);    % constant and exogenous variables (hd)
    HDshocks=zeros(info.nvar,info.horizons,info.nvar,info.nfinal); % shocks (hd)
    HDbaseline=zeros(info.nvar,info.horizons,info.nfinal);         % together: lag-structure + constant +exogenous variables (hd)

    for i=1: info.nfinal
        draw = randsample(1 :count_accept, 1,true,nisw(1:count_accept,1));
        BDraw= BDraws{draw,1};
        sigmaDraw= SigmaDraws{draw,1};
        QDraw= QDraws{draw,1} ;
        x  = [vec(BDraw); vec(sigmaDraw);vec(QDraw)];
        y  =f_h_map_inv(x,info);
        % calculate all impulse responses
        L=irf(y,info); 
            irf_matrix=zeros(info.nvar,info.nvar,info.horizons);
            for j=0:info.horizons-1
                irf_matrix(:,:,j+1)=L(int2str(j));
            end
         % calculate hystorical decomposition (hd)
        [hd_a,hd_b,hd_c]=hd(y,irf_matrix,info);    
        IrfDraws(:,:,:,i)=irf_matrix;
        HDaDraws(:,:,i)=hd_a;
        HDbDraws(:,:,:,i)=hd_b;
        HDshocks(:,:,:,i)=hd_c; 
    end

    %add lag-structure, constant and exogenous variables to baseline contribution
    HDbaseline=HDaDraws+squeeze(sum(HDbDraws,3));


    %*****************************************************************************************************************
    %% PLOT IMPULSE RESPONSES AND HYSTORICAL DECOMPOSITION
    %*****************************************************************************************************************

    % Use {plotIrf} function to plot Impulse Response Functions

    % INPUT: 
    % irfArray:       calculated IRFs (variable, shock, horizont, draw)
    % info:           VAR info
    % shockList:      list of shocks, for example: 1:4 or [1,3]
    % variableList:   list of variables, for example: 1:6, [1,4,5]
    % cumulativeList: list  of variables for which the corresponding impulse responses need to be accumulated

    figure('units','normalized','position',[.1 .1 .95 .95])
    plotIrf(IrfDraws,info,[1,2,3,4],[1,2,3,4],[1,2,3,4]);
    irf_name=strcat(sheet,"_irf.pdf");
    save2pdf(irf_name);
    close;


    %**************************************************************************************************************

    % Use {plotHd} function to plot Historical Decomposition

    % INPUT:
    % HDbaseline:      (variable,horizon,shock) /includes: lag-structure, constant and exogenous variables contribution
    % HDshocks:        (variable,horizon,shock,draw) / shocks contribution
    % info:            VAR information
    % varSalection:    select variable for drawing historical decomposition
    % shockSelection:  select shocks to plot
    % accumP:          select number of periods to accumulate
    %                  (To have YoY growth rates when using log(diff) transformation accumP is 4)
    % addBaseline:     add baseline contribution to figure('yes' or 'no')
    % colors:          set the colors for historical decomposition bar plot
    %                  last color is for baseline and then other shocks { ...,[0 0 1],..}
    %                  red ([1 0 0]), green([0 1 0]), blue([0 0 1]), yellow([1 1 0]), black([0 0 0 ])

    colors={[1 0 0],[0 1 0],[0 0 1],[0.5 0.5 0.5],[1 1 0]};
    addBaseline='yes';
    accumP=[4,4,4,4];
    shockSelection=[1,2,3,4];
    varSelection=[1,2,3,4];

    for i=varSelection
         figure('units','normalized','position',[.1 .1 .95 .95])
          plotHd(HDbaseline,HDshocks,info,i,shockSelection,accumP(i),addBaseline,colors);
          hd_name=strcat(sheet,'_',info.names{i},'_hd.pdf');
          save2pdf(hd_name);
          close;
    end

%*******************************************************************************************************
% SAVE  
%*******************************************************************************************************

%save workspace
workspace_name=strcat(sheet,'_workspace.mat');
save(workspace_name);

% save output of "draw_nwp" or "draw_inwp" function to excel file  
excel_name=strcat(sheet,'_output.xlsx');
saveOutputExcel(excel_name,BDraws,SigmaDraws,QDraws,WDraws,nisw,uisw,vol1,vol2,count_accept,nfinal)

% save irfs to excel file (lower bound, median and upper bound)
excel_name=strcat(sheet,'_output.xlsx');
saveIrfExcel(excel_name,info,IrfDraws)
%*********************************************************************************************************
    
    
end


