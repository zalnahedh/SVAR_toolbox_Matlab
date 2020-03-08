function test = testSpd(X)
% function tests if matrix X is symmetric positive definite(SPD)
% output: 1 -> matrix X is SPD
%         0 -> matrix X is not SPD
 
% When flag output is specified, "chol" function does not generate an error if the input matrix is not symmetric positive definite.
% If flag = 0 then the input matrix is SPD and the factorization was successful
% If flag flag >=1 , then the input matrix is not SPD and flag is an integer indicating the index of the pivot position where the factorization failed
[~,flag]=chol(X,'lower');

if flag==0 
    test=1 
else
    test=0
end

end

