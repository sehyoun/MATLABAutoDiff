function x = rdivide(x,y)
% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at

if isa(x, 'myAD')
    if isa(y, 'myAD')
        x.derivatives = valXder(1./y.values, x.derivatives) - valXder(x.values./y.values.^2, y.derivatives);
        x.values = x.values./y.values;
    else
        x.derivatives = valXder(1./y, x.derivatives);
        x.values = x.values./y;
    end
else
    y.derivatives = valXder(- x./y.values.^2, y.derivatives);
    y.values = x./y.values;
    x = y;
end
