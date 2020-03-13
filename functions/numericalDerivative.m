function Df = numericalDerivative(f, x, epsilon,order)
%-------------------------------------------------------------------------
% INPUT:  f        - function from R^n into R^m
%         x        - a vector of length n
%         epsilon  - small change in x
%         order    - 2nd or 4th order approximation
% OUTPUT: Df       - m x n  matrix of partial derivatives evaluated at x
%                    Df(i,j) = df_i/dx_j
%------------------------------------------------------------------------

n=numel(x);
% if epsilon is missing put some default value
if nargin < 3
    epsilon=1.0e-6;
    order=2;
end
z=x;

if (order==2)
    % option 1 (2nd order approximation)
    for j=1:n
        z(j)=x(j)+epsilon;
        y1=f(z);
        z(j)=x(j)-epsilon;
        y0=f(z);
        z(j)=x(j);
        if j == 1
            Df=zeros(numel(y0),n);
        end
        Df(:,j)=(y1-y0)/(2*epsilon);
    end
end

% option 2 (4th order approximation)

if (order==4)
    for j=1:n
        z(j)=x(j)+epsilon;
        yp1=f(z);
        z(j)=x(j)+2*epsilon;
        yp2=f(z);
        z(j)=x(j)-epsilon;
        ym1=f(z);
        z(j)=x(j)-2*epsilon;
        ym2=f(z);
        z(j)=x(j);
          if j == 1
            Df=zeros(numel(ym1),n);
          end
        Df(:,j)=(-yp2+8*yp1-8*ym1+ym2)/(12*epsilon);
    end
end

end