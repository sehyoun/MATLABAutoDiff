function x = minus(x,y)
if isa(x, 'myAD')
    if isa(y, 'myAD')
        x.values = x.values - y.values;
        x.derivatives = x.derivatives - y.derivatives;
    else
        x.values = x.values - y;
    end
else
    y.values = x - y.values ;
    y.derivatives = - y.derivatives;
    x = y;
end
