function varargout = sort(x)
% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at

[val, idx] = sort(x.values);
x.values = val;
x.derivatives = x.derivatives(idx,:);
varargout{1} = x;
if (nargout>1)
    varargout{2} = idx;
end
