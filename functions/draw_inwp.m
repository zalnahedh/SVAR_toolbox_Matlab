function [BDraws,SigmaDraws,QDraws,nisw,count_accept,nfinal] = draw_inwp(info,d2)


    count_accept_guess      = 1000;
    BDraws                  = cell([count_accept_guess ,1]); % reduced-form  parameters
    SigmaDraws              = cell([count_accept_guess ,1]); % reduced-form covariance matrices
    QDraws                  = cell([count_accept_guess ,1]); % orthogonal matrices
    v1_volume_element       = zeros(count_accept_guess ,1);  % absolute value of logarithm of overall volume element in numerator
    v2_volume_element       = zeros(count_accept_guess ,1);  % absolute value of logarithm of overall volume element in denumerator
    uisw                    = zeros(count_accept_guess ,1);  % unnormalized importance sampler weights


    % start Gibbs sampling from OLS estimate of beta
    beta=vec(invpd((info.X)'*info.X)*(info.X)'*(info.Y));
    
    
    % counters
    j_burn=1;             %count number of rejected draws until convergence of Markov chain 
    record=1;             %count number of importance sampling draws
    count_iterations = 1; %count number of importance sampling draws before showing current status
    count_accept   = 0;   %count number of draws satisfying zero/sign restrictions
    nfinal=0;             %initialize effective sample size
    while ((record<=info.maxDraws) && (nfinal<=info.finalDraws)) 
            %display(['DRAW = ',num2str(record)])
            %tic
            sigmaDraw=drawSigmaCondBeta(info,beta,d2.scaleMatrix0,d2.alfa0);
            stability=0;
            % checking SVAR stablity 
            while (stability==0)
                BDraw=drawBetaCondSigma(info,d2.beta0,d2.omega0,sigmaDraw);
                stability=checkVarStability(info,BDraw);
            end   
            beta=BDraw;
            BDraw = reshape(BDraw,info.m,info.nvar);
            j_burn=j_burn+1;    
            
            if j_burn>d2.burn
                     % to speed up proccess draw several Q matrices for the same beta and sigma draws
                     for jj=1:info.maxQ
                            %display('  Drawing vectors w...');
                            w  = drawVectorsW(info);
                            x  = [vec(BDraw); vec(sigmaDraw);w];

                            % gf_h_map_inv is mapping from (Beta,Sigma,W) to structural representation (A0,Aplus)
                            % by construction this structural representation satisfies  zero restrictions  
                            % display('  Calculating structural representation that satisfies zero restrictions...');
                            y  =gf_h_map_inv(x,info);

                            % Q is already calculated in gf_h_map_inv to obtain final structural representation that satisfies zero restrictions
                            %if we want to save Q matrix wToQ function is used
                            %display('  Calculating Q matrix...');
                            QDraw = wToQ(x,info);

                            % store reduced-form draws and Q matrix (including not accepted draws)
                            %BDraws{record,1}          = BDraw;
                            %SigmaDraws  {record,1}    = sigmaDraw;
                            %QDraws{record,1}          = reshape(QDraw,nvar,nvar);

                            % calculate unnormalized importance sampler weights(uisw) if structural representation satisfies sign restrictions
                            % otherwise set the unnormalized importance sampler weights to zero
                            if checkSignRes(y,info)
                                    %display('  Sign restrictions satisfied, calculating importance sampler weights...');
                                    count_accept=count_accept+1;
                                    BDraws{count_accept,1}        = BDraw;
                                    SigmaDraws  {count_accept,1}  = sigmaDraw;
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
                                    display(['Remaining number of reduced_form parameters draws = ',num2str(info.maxDraws-(record))])
                                    display(['Number of accepted parameters = ',num2str(count_accept)])
                                    display(['Number of nfinal parameters = ',num2str(nfinal)])
                                    count_iterations =0;
                           end
                           count_iterations = count_iterations + 1;
                           record=record+1;
                     end
            end
    end
BDraws             =  BDraws(~cellfun('isempty',BDraws ));            
SigmaDraws         =  SigmaDraws(~cellfun('isempty',SigmaDraws)); 
QDraws             =  QDraws(~cellfun('isempty',QDraws ));  
nisw               =  nisw(1:count_accept);
    
end

