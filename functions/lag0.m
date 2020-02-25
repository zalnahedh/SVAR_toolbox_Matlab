function output=lag0(x,pp)
 
%---Desription--------------------------------------------------------
%For an input vector or matrix x=[a_{1}; a_{2}; ...; a_{n}] where a_{i}%
%are row vectors, it returns the vector or matrix:                     %
%                                                                      % 
%output= [0; ... ; 0 ; a_{1}; a_{2}; ....; a_{n-p}]                      %
%                                                                      % 
%In other words it lags the variable p periods and places zeros in the % 
%rows corresponding to the first p periods                             % 
%                                                                      %                                                     % 
%----------------------------------------------------------------------


%Compute the number of rows and columns of input x and save in 
%variables R and C respectively 
[R,C]=size(x);
 
%Take the first R-p rows of matrix x and all columns
x1=x(1:(R-pp),:);
%Preceed them with p rows of zeros and return
output=[zeros(pp,C); x1];
 
end

