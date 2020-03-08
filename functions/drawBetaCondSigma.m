function [beta] = drawBetaCondSigma(info,beta0,omega0,sigma)

q=info.nvar*(info.nvar*info.nlag+info.nex+info.cte);

%omega0Inv=invpd(omega0);

% omega0 is diagonal matrix, inverse is easily calculated
omega0Inv=diag(1./diag(omega0));

sigmaInv=invpd(sigma);
omega1=invpd(omega0Inv+kron(sigmaInv,(info.X)'*info.X));
beta1=omega1*(omega0Inv*beta0+kron(sigmaInv,(info.X)')*vec(info.Y)); 
omega1=(omega1+omega1')/2; %symmetric
beta=beta1+(randn(1,q)*chol(omega1))';

end

