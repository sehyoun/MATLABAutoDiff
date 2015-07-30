function z = min(varargin)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

    if nargout > 1
        error('Only one output for min implemented');
    end

    if nargin > 1
        x = min(varargin{1}.values, varargin{2}.values);
        idx = (x == varargin{2}.values) + 1;
        z = varargin{idx};
    else
        [x, idx] = min(varargin{1}.values);
        z = varargin{1};
        z.values = z.values(idx);
        z.derivatives = z.derivatives(idx,:);
    end
