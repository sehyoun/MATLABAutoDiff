function x = times(x,y)
if isa(x, 'myAD')
    l=size(x.derivatives,2);
    if isa(y, 'myAD') 
        x.derivatives = repmat(y.values(:),1,l).*x.derivatives + repmat(x.values(:),1,l).*y.derivatives;
        x.values = x.values.*y.values;
    else
        x.values = x.values.*y;
        x.derivatives = repmat(y(:),1,l).*x.derivatives;
    end
else
    y.values = x.*y.values;
    y.derivatives = repmat(x(:),1,size(y.derivatives,2)).*y.derivatives;
    x = y;
end