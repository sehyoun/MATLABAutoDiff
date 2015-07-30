function x=transpose(x)
    [n,m]=size(x.values);
    p=reshape(reshape(1:n*m,m,n)',n*m,1);
    x.derivatives=x.derivatives(p,:);
    x.values=x.values.';