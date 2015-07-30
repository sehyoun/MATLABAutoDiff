function y = subsasgn(y, S, x)
%% No Subsetting %%
if  isempty(S.subs{1})
    return;
end

if isa(x, 'myAD')
    if isa(y, 'myAD')
        y.values(S.subs{:}) = x.values;
        [n,m]=size(y.values);
        locs=reshape(1:n*m,n,m);
        locs=locs(S.subs{:});
        y.derivatives(locs(:),:) = x.derivatives;
    else
        [n,m]=size(y);
        locs=reshape(1:n*m,n,m);
        locs=locs(S.subs{:});
        l = size(x.derivatives,2);
        y = myAD(y,sparse(n*m,l));
        y.values(S.subs{:}) = x.values;
        y.derivatives(locs(:),:) = x.derivatives;
    end
else
    [n,m]=size(y.values);
    locs=reshape(1:n*m,n,m);
    locs=locs(S.subs{:});
    y.values(S.subs{:}) = x;
    y.derivatives(locs(:),:) = 0;
end