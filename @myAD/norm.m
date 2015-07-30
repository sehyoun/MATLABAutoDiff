function x = norm(x, p)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

    if (nargin==1) p = 2; end
    x = sum(abs(x).^p).^(1/p);
