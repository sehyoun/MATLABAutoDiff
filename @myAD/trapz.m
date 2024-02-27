function Q = trapz(varargin)
  if nargin == 1
    Y = varargin{1};
    n = length(Y);  
    dx = ones(1, n-1);
  elseif nargin == 2
    Y = varargin{2};
    n = length(Y);  
    dx = diff(varargin{1});
    dx = reshape(dx, 1, n-1);
  else
    error("not implemented yet")
  end
  A = spdiags(ones(n-1, 1), 0, n-1, n) + spdiags(ones(n-1,1), 1, n-1, n);
  dxA = dx*A;
  Y = reshape(Y, n, 1);
  Q = dxA*Y/2;
end
