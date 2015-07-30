function x = sum(x,varargin)
if nargin<2
    [n,m]=size(x.values);
    l=size(x.derivatives,2);
    x.derivatives = reshape(sum(reshape(x.derivatives,n,m*l)),m,l);
    x.values = sum(x.values);
elseif varargin{1}==2
    x=sum(x');
end