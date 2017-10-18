function x = minus(x,y)
% Updated by SeHyoun Ahn, Oct 2017

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if isa(x, 'myAD')
    if isa(y, 'myAD')
        x.values = x.values - y.values;
        if numel(x) == numel(y)
            x.derivatives = x.derivatives - y.derivatives;
        elseif (numel(y) == 1) || (numel(x) == 1)
            x.derivatives = bsxfun(@minus, x.derivatives, y.derivatives);
        else
            error('Implicit expansion is not supported yet!');
        end
    else
        x.values = x.values - y;
    end
else
    y.values = x - y.values ;
    y.derivatives = - y.derivatives;
    x = y;
end
