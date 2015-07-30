function idx = end(x, k, n)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

if (k==1)
    idx = length(x.values);
else
    idx = 1;
end
