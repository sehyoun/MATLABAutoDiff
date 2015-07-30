function x = cumprod(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

    n = size(x.values,1);
    for i = 2:n
        x.derivatives(i,:) = x.values(i-1).*x.derivatives(i,:) + x.values(i).*x.derivatives(i-1,:);
        x.values(i) = x.values(i-1).*x.values(i);
    end
