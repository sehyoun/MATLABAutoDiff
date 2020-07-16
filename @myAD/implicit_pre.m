function x0 = implicit_pre(x0, l_pre)
  % by SeHyoun Ahn, April 2020

  % Get number of derivatives
  if nargin < 2
    l_pre = size(getderivs(x0), 2);
  elseif isa(l_pre, 'myAD')
    l_pre = size(getderivs(l_pre), 2);
  end

  % getvalues if AD
  if isa(x0, 'myAD')
    x0 = getvalues(x0);
  end

  n_x = length(x0);
  x0 = myAD(x0, spdiags(ones(n_x,1), l_pre, n_x, n_x+l_pre));
end
