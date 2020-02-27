function [beta0,omega0] = minnesotaPrior(info,lambda,blockExo)

% returns prior values from hyperparameters, for the Minnesota prior
% inputs:  
%          - scalar 'lambda1': overall tightness hyperparameter (defined p 16 of technical guide)
%          - scalar 'lambda2': cross-variable weighting hyperparameter(defined p 16 of technical guide)
%          - scalar 'lambda3': lag decay hyperparameter (defined p 16 of technical guide)
%          - scalar 'lambda4': exogenous variable tightness hyperparameter (defined p 17 of technical guide)
%          - scalar 'lambda5': block exogeneity shrinkage hyperparameter (defined p 32 of technical guide)
%          - integer 'n': number of endogenous variables in the BVAR model (defined p 7 of technical guide)
%          - integer 'm': number of exogenous variables in the BVAR model (defined p 7 of technical guide)
%          - integer 'p': number of lags included in the model (defined p 7 of technical guide)
%          - integer 'k': number of coefficients to estimate for each equation in the BVAR model (defined p 7 of technical guide)
%          - integer 'prior': value to determine which prior applies to the model
%          - integer 'bex': 0-1 value to determine if block exogeneity is applied to the model
%          - matrix 'blockexo': matrix indicating the variables on which block exogeneity must be applied
% outputs: - vector 'beta0': vector of prior values for beta (defined in 1.3.4) 
%          - matrix 'omega0': prior covariance matrix for the VAR coefficients (defined in 1.3.8)
%          - matrix 'sigma': 'true' variance-covariance matrix of VAR residuals, for the original Minnesota prior


% VAR characteristics:

m=info.nex+info.cte;
k=info.nvar*info.nlag+m; %total number of VAR coefficients per each equation
q=info.nvar*k; %total number of VAR coefficients

% Minnesota prior hyperparameters
L1=lambda.lambda1;
L2=lambda.lambda2;
L3=lambda.lambda3;
L4=lambda.lambda4;
L5=lambda.lambda5;

firstLagCoef=lambda.firstLagCoef;


%block exogeneity
bex=blockExo.bex;
bexMatrix=blockExo.bexMatrix;

%initialization
beta0=zeros(q,1);
omega0=zeros(q,q);


% initialize beta0 (mean):

beta0(1,1)=firstLagCoef;
for ii=1:info.nvar
    beta0((ii-1)*k+ii,1)=firstLagCoef;
end


% set the variance on coefficients trelated to own lags
for ii=1:info.nvar
   for jj=1:info.nlag
   omega0((ii-1)*k+(jj-1)*info.nvar+ii,(ii-1)*k+(jj-1)*info.nvar+ii)=(L1/jj^L3)^2;
   end
end

%  set variance for coefficients on cross lags
for ii=1:info.nvar
   for jj=1:info.nlag
      for kk=1:info.nvar
      if kk==ii
      else
      omega0((ii-1)*k+(jj-1)*info.nvar+kk,(ii-1)*k+(jj-1)*info.nvar+kk)=((lambda.stand(ii)/lambda.stand(kk))^2)*(((L1*L2)/(jj^L3))^2);
      end
      end
   end
end

% set the variance for exogenous variables, using (1.3.7)
for ii=1:info.nvar
   for jj=1:m
   omega0(ii*k-m+jj,ii*k-m+jj)=(lambda.stand(ii)*L1*L4)^2;
   end
end

% if block exogeneity has been selected, implement it, according to (1.7.4)
if bex==1
   for ii=1:info.nvar
      for jj=1:info.nvar
         if bexMatrix(jj,ii)==1
            for kk=1:info.nlag
            omega0((jj-1)*k+(kk-1)*info.nvar+ii,(jj-1)*k+(kk-1)*info.nvar+ii)=omega0((jj-1)*k+(kk-1)*info.nvar+ii,(jj-1)*k+(kk-1)*info.nvar+ii)*L5^2;
            end
         else
         end
      end
   end
% if block exogeneity has not been selected, don't do anything 
else
end

end

