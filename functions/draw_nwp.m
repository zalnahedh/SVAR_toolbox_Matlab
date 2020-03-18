function [BDraws,SigmaDraws,QDraws,WDraws,nisw,uisw,vol1,vol2,count_accept,nfinal] = draw_nwp(info,d1)
%******************************************************************************************************
%   INPUT:
%       info ->  structure / SVAR information
%       d1   ->  structure / impose prior (NORMAL WISHART PRIOR) 

%   OUTPUT:
%       BDraws          ->  all reduced-form  parameters
%       SigmaDraws      ->  all reduced-form covariance matrices
%       QDraws          ->  all orthogonal matrices
%       WDraws          ->  all w draws
%       nisw            ->  normalized importance sampling weights
%       uisw            ->  unnormalized importance sampling weights
%       vol1            ->  absolute value of logarithm of overall volume element in numerator
%       vol2            ->  absolute value of logarithm of overall volume element in denumerator
%       count_accept    ->  total number of draws satisfying restrictions
%       nfinal          ->  effective number of draws satisfying restrictions
%**************************************************************************************************



%=====================================================================================================================================================
%% INITIALIZATION

% we don't know effective size in advance (initial guess for count_accept needed)
count_accept_guess      = 20000;

BDraws                  = cell([count_accept_guess ,1]); % reduced-form  parameters
SigmaDraws              = cell([count_accept_guess ,1]); % reduced-form covariance matrices
QDraws                  = cell([count_accept_guess ,1]); % orthogonal matrices
WDraws                  = cell([count_accept_guess ,1]); % w draws
v1_volume_element       = zeros(count_accept_guess ,1);  % absolute value of logarithm of overall volume element in numerator
v2_volume_element       = zeros(count_accept_guess ,1);  % absolute value of logarithm of overall volume element in denumerator
uisw                    = zeros(count_accept_guess ,1);  % unnormalized importance sampler weights

    record=1;               %count number of importance sampling draws
    count_iterations = 1;   %count number of importance sampling draws until showing current status
    count_accept   = 0;     %count number of draws satisfying restrictions
    nfinal=0;               %initialize effective sample size
    
    wb=waitbar(0);
    
    while ((record<=info.maxDraws) && (nfinal<=info.finalDraws))
            %display(['DRAW = ',num2str(record)])
            %tic
            %draw sigma from IW-distribution parametrized with scale matrix phiTilda and nuTilda degrees of freedom
            sigmaDraw = iwpq(d1.nuTilda,inv(d1.phiTilda));
            %make Cholesky factorization
            cholSigmaDraw = (info.h(sigmaDraw))';
            
            stability=0;
            % checking stablity influence distributions?????
            while stability==0    
                BDraw = kron(cholSigmaDraw,d1.cholOmegaTilda)*randn(info.m*info.nvar,1) + reshape(d1.psiTilda,info.nvar*info.m,1);
                stability=checkVarStability(info,BDraw);
            end   
            BDraw = reshape(BDraw,info.m,info.nvar);
 
            for jj=1:info.maxQ
                    %display('  Drawing vectors w...');
                    w  = drawVectorsW(info);
                    x  = [vec(BDraw); vec(sigmaDraw);w];

                    % gf_h_map_inv is mapping from (Beta,Sigma,W) to structural representation (A0,Aplus)
                    % and by construction this structural representation satisfies  zero restriction:  
                    %display('  Calculating structural representation that satisfies zero restrictions...');
                    y  =gf_h_map_inv(x,info);

                    % Q is already calculated in gf_h_map_inv to obtain final structural representation that satisfies zero restriction, but if we want to
                    % save Q matrix wToQ function is used
                    %display('  Calculating Q matrix...');
                    QDraw = wToQ(x,info);

                    % store reduced-form draws and Q matrix (including not accepted draws)
                    %BDraws{record,1}          = BDraw;
                    %SigmaDraws  {record,1}    = sigmaDraw;
                    %QDraws{record,1}          = reshape(QDraw,nvar,nvar);

                    % calculate unnormalized importance sampler weights(uisw) if structural representation satisfies sign restrictions
                    % otherwise set the unnormalized importance sampler weight to zero
                    if checkSignRes(y,info)
                            %display('  Sign restrictions satisfied, calculating importance sampler weights...');
                            count_accept=count_accept+1;
                            BDraws{count_accept,1}        = BDraw;
                            SigmaDraws  {count_accept,1}  = sigmaDraw;
                            WDraws{count_accept,1}        = w;
                            QDraws{count_accept,1}        = reshape(QDraw,info.nvar,info.nvar);
                            [v1,v2,weight] = computeWeight(y,info,info.type_of_restrictions);
                            %toc
                            v1_volume_element(count_accept,1) = v1;
                            v2_volume_element(count_accept,1) = v2;
                            uisw(count_accept,1)              = weight; 
                            %effective size
                            nisw=uisw/sum(uisw);
                            ne=floor(1/sum(nisw.^2));
                            nfinal=min(ne,count_accept);
                            if (nfinal>=info.finalDraws) 
                                display(['Break.Number of accepted parameters = ',num2str(count_accept)])
                                display(['Break.Number of nfinal parameters = ',num2str(nfinal)])
                                break;
                            end
                   else
                            %display('  Sign restrictions NOT satisfied!...');
                            %toc
                            %v1_volume_element(record,1) = NaN;
                            %v2_volume_element(record,1) = NaN;
                            %uisw(record,1) = 0.00000;     
                   end

                   if count_iterations == info.iter_show
                            display(['Number of reduced-form parameters draws = ',num2str(record)])
                            display(['Maximum number of remaining number of reduced_form parameters draws = ',num2str(info.maxDraws-(record))])
                            display(['Number of accepted parameters = ',num2str(count_accept)])
                            display(['Number of nfinal parameters = ',num2str(nfinal)])
                            count_iterations =0;
                   end
                   waitbar(nfinal/info.finalDraws, wb, sprintf('    Effective sample: %d  ; Number of accepted parameters: %d     ',nfinal, count_accept));
                   count_iterations = count_iterations + 1;
                   record=record+1;
            end
    end
 
% output     
BDraws             =  BDraws(~cellfun('isempty',BDraws ));            
SigmaDraws         =  SigmaDraws(~cellfun('isempty',SigmaDraws)); 
QDraws             =  QDraws(~cellfun('isempty',QDraws ));
WDraws             =  WDraws(~cellfun('isempty',WDraws ));
nisw               =  nisw(1:count_accept);
uisw               =  uisw(1:count_accept);
vol1               =  v1_volume_element(1:count_accept);
vol2               =  v2_volume_element(1:count_accept);
    
    
end

