function x = mtimes(x, y)

%% MATRIX MULTIPLICATION %%
if isa(x,'myAD')
    n=size(x.values,1);
    m=size(x.values,2);
    l=size(x.derivatives,2);
    if isa(y,'myAD')
        if m>1 && length(y)==m
            x.derivatives=x.derivatives.*repmat(y.values,n,l)+repmat(y.derivatives,n,1).*repmat(reshape(x.values',n*m,1),1,l);
            x.derivatives=reshape(sum(reshape(x.derivatives,n,m*l)),n,l);
            x.values=x.values*y.values;
        end
    else
        x.derivatives=x.derivatives.*repmat(y,n,l);
        x.derivatives=reshape(sum(reshape(x.derivatives,n,m*l)),n,l);
        x.values=x.values*y;        
    end
else
    n=size(x,1);
    m=size(x,2);
    l=size(y.derivatives,2);
    y.derivatives=repmat(y.derivatives,n,1).*repmat(reshape(x',n*m,1),1,l);
    y.derivatives=reshape(sum(reshape(y.derivatives,n,m*l)),n,l);
    y.values=x*y.values;
    x=y;
end

return;
