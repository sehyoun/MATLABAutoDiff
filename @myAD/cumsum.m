function x = cumsum(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

    x.derivatives = cumsum(x.derivatives);
    x.values = cumsum(x.values);
