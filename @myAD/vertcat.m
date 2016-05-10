function x = vertcat(varargin)
% Edited by SeHyoun Ahn, Jan 2016
% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at
i = 1;
while (~isa(varargin{i}, 'myAD'))
    i=i+1;
end
y=vertcat(varargin{1:i-1});
n=size(y,1);
x = varargin{i};
l=size(x.derivatives,2);
[nx,m]=size(x.values);
if (i>1)
    x.values = [y; x.values];
    x.derivatives = [sparse(n,m*l); reshape(x.derivatives,nx,m*l)];
else
    x.derivatives = reshape(x.derivatives,nx,m*l);
end
n=nx+n;

for j = i+1:nargin
    if isa(varargin{j}, 'myAD')
        x.values = [x.values; varargin{j}.values];
        nvar=size(varargin{j}.values,1);
        n=n+nvar;
        x.derivatives = [x.derivatives; reshape(varargin{j}.derivatives,nvar,m*l)];
    else
	n=n+size(varargin{j},1);
        x.values = [x.values; varargin{j}];
        x.derivatives(end+size(varargin{j},1),:) = 0;
    end
end
x.derivatives=reshape(x.derivatives,n*m,l);
