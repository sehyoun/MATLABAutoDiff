function x=ctranspose(x)
    [n,m]=size(x.values);
    p=reshape(reshape(1:n*m,n,m)',n*m,1);
    x.derivatives=x.derivatives(p,:);
    x.values=x.values';