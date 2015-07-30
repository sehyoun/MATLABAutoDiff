function x = vertcat(varargin)
y = []; i = 1;
while (~isa(varargin{i}, 'myAD'))
    y = [y; varargin{i}];
    i=i+1;
end
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
        x.values = [x.values, varargin{j}];
        x.derivatives(end+size(varargin{j},1),:) = 0;
    end
end
x.derivatives=reshape(x.derivatives,n*m,l);