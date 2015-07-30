function x = diff(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

    x.values = diff(x.values);
    x.derivatives = diff(x.derivatives);
