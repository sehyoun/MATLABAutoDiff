function display(x, varargin)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

if (nargin < 2)
    disp([inputname(1) ':']);
else
    disp([varargin{1} ':']);
end
disp('Values =');
disp(x.values);
disp('Derivatives =');
disp(x.derivatives);
