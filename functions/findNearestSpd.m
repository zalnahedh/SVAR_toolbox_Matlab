function A_new = findNearestSpd(A)
% function computes the nearest SPD matrix to A,
% by methodology proposed by N.J Highams (1988): 'Computing a nearest symmetric positive semidefinite matrix',Linear Algebra and its Applications, No. 103

%check whether matrix is SPD or not
test1=testSpd(A);

if (test1==1)
    A_new=A;
else
    B=(A+A')/2;
    [~,S,V]=svd(B);
    H=V*S*V';
    new=(B+H)/2;
    new=(new+new')/2;
    test2=testSpd(new);
    if (test2==1)
        A_new=new;
    else  
       n=1;
       while (test2==0) %small modifications to achieve SPD 
          minEig=min(eig(new));
          new=new+(-minEig*n.^2+eps(minEig)+0.0000000000000000001)*eye(size(new));  
          test2=testSpd(new);
          n=n+1;
       end
       A_new=new;   
     end
end



end

