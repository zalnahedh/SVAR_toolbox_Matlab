function [v1,v2,weight] = computeWeight(x,info,s)
%=================================================================================================================
nvar   =info.nvar;
m      =info.m;
A0     =reshape(x(1:nvar*nvar),nvar,nvar);

strToW            =info.strToW;
zeroRestr         =info.zeroRestr;
strToRed          =info.strToRed;
strToIrf          =info.strToIrf;
irfToStr          =info.irfToStr;


 switch s
     case 'structural'
            v1      =   (nvar*(nvar+1)/2)*log(2)-(2*nvar+m+1)*logAbsDet(A0);
            v2      =   logVolumeElement(strToW,x,zeroRestr); 
            weight  =   exp(v1-v2);

     case 'irf'
           % Transform 
           irf_param    =   strToIrf(x);
           % we first transform parameters from reduced-form to structural and then from structural to irf representation
           % so, we need to calculate product of two volume elements (sum of log volume elements)
           % First: volume element of structural to reduced-form mapping calculated at structural representation point x
           % Second: volume element of irf to structural mapping calculated at irf representation point irf_param
           v1           =   logVolumeElement(strToRed,x) + logVolumeElement(irfToStr, irf_param);
           % For proposal density calculation of volume elements is little more complicated and logVolumeElement function takes 3 arguments
           % Third argument is zeroRestrictions function
           v2           =   logVolumeElement(strToW,x,zeroRestr) + logVolumeElement(irfToStr, irf_param ,zeroRestr); 
           weight       =   exp(v1-v2);
    
     case 'both'
         
         return 
                         
     otherwise
             v1     =1;
             v2     =1;
             weight =1;
             display(' Weight is set to 1. Please select one of strings (irf, structural or both) as third argument.');
 end

end

