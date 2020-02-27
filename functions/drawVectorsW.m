function w = drawVectorsW(info)
 
nvar=info.nvar;
dimSphere=info.dimSphere;
W_matrix=info.W_matrix;
 
w=zeros(dimSphere,1);
k=1;
for j=1:nvar
    s=size(W_matrix{j},1);
    wj=randn(s,1);
    w(k:k+s-1)=wj/norm(wj);
    k=k+s;
end
