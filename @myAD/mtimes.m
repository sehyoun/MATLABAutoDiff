function x = mtimes(x, y)

%% MATRIX MULTIPLICATION %%
if isa(x,'myAD')
    [n,m]=size(x.values);
    l=size(x.derivatives,2);
    if isa(y,'myAD')
        if m>1 && length(y)==m
            locs=reshape(1:n*m,n,m)';
            x.derivatives=x.derivatives(locs,:);
            x.derivatives=x.derivatives.*repmat(y.values,n,l)+repmat(y.derivatives,n,1).*repmat(reshape(x.values',n*m,1),1,l);
            x.derivatives=reshape(sum(reshape(x.derivatives,n,m*l)),n,l);
            x.values=x.values*y.values;
        elseif max(m,n)==1
            x.derivatives=repmat(x.derivatives,numel(y.values),1).*repmat(y.values(:),1,l)+y.derivatives*x.values;
            x.values=x.values*y.values;
        elseif numel(y.values)==1
            x.derivatives=x.derivatives*y.values+repmat(x.values,1,l).*repmat(y.derivatives,n*m,1);
            x.values=x.values*y.values;
        else
            error('This form of multiplication not supported');
        end
    else
        if m>1 && length(y)==m
            locs=reshape(1:n*m,n,m)';
            x.derivatives=x.derivatives(locs,:)';
            x.derivatives=x.derivatives.*repmat(y,n,l);
            x.derivatives=reshape(sum(reshape(x.derivatives,n,m*l)),n,l);
            x.values=x.values*y;
        elseif max(m,n)==1
            x.derivatives=repmat(x.derivatives,numel(y),1).*repmat(y(:),1,l);
            x.values=x.values*y;
        elseif numel(y)==1
            x.derivatives=x.derivatives*y;
            x.values=x.values*y;
        else
            error('This form of multiplication not supported');
        end
    end
else
    n=size(x,1);
    m=size(x,2);
    l=size(y.derivatives,2);
    if m>1 && length(y)==m
        y.derivatives=repmat(y.derivatives,n,1).*repmat(reshape(x',n*m,1),1,l);
        y.derivatives=reshape(sum(reshape(y.derivatives,n,m*l)),n,l);
        y.values=x*y.values;
        x=y;
    elseif max(m,n)==1
        y.derivatives=y.derivatives*x;
        y.values=x*y.values;
        x=y;
    elseif numel(y.values)==1
        y.derivatives=repmat(x,1,l).*repmat(y.derivatives,n*m,1);
        y.values=x*y.values;
        x=y;
    else
        error('This form of multiplication not supported');
    end
end