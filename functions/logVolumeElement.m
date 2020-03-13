function vol_elem = logVolumeElement(f,x,h)
%
%  Returns the log of the volume element of f restricted to the set of points x
%  in R^n that satisfy h(x) = 0.
%
%  f is a function from R^n to R^m.
%  x is a point in R^n satisfying h(x)=0.
%  h is a function from R^n to R^(n-k).
%  
%  It is assumed that Dh(x) is of full row rank whenever h(x) = 0.  Under this
%  assumption h defines a k-dimensional manifold in R^n.  ve is the volume
%  element of f, restricted to the manifold, when evaluated at x.
%

% change numerical derivative parameters here
epsilon=1.0e-6;
order=2;

Dfx=numericalDerivative(f,x,epsilon,order);          % m x n matrix
if nargin > 2
    Dhx=numericalDerivative(h,x,epsilon,order);      % (n-k) x n matrix
    N=null(Dhx);
    M=Dfx*N;                  % perp(Dhx') - n x k matrix
    vol_elem=0.5*logAbsDet(M'*M);
else
    vol_elem=0.5*logAbsDet(Dfx'*Dfx);
end
end
