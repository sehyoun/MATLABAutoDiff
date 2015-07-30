function x = abs(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

    tmp = sign(x.values);
    x.derivatives = tmp(:,ones(size(x.derivatives,2),1)).*x.derivatives;
    x.values = abs(x.values);