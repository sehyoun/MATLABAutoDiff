function x = plus(x, y)
% Updated by SeHyoun Ahn, Oct 2017

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if isa(x, 'myAD')
    if isa(y, 'myAD')
        x.values = x.values + y.values;
        if numel(x) == numel(y)
            x.derivatives = x.derivatives + y.derivatives;
        elseif (numel(x) == 1) || (numel(y) == 1)
            x.derivatives = bsxfun(@plus, x.derivatives, y.derivatives);
        else
            error('Implicit Expansion is not supported yet!');
        end
    else
        x.values = x.values + y;
    end
else
    y.values = x + y.values;
    x = y;
end
