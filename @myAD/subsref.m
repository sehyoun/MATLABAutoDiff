function x = subsref(x, S)
    [n,m]=size(x.values);
    locs=reshape(1:n*m,n,m);
    locs=locs(S.subs{:});
    x.derivatives = x.derivatives(locs(:),:);
    x.values = x.values(S.subs{:});