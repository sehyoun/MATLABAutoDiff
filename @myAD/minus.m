function x = minus(x,y)
  % Updated by SeHyoun Ahn, Oct 2017

  % In Package myAD - Automatic Differentiation
  % by Martin Fink, June 2006
  % martinfink 'at' gmx.at
  if isa(x, 'myAD')
    if isa(y, 'myAD')
      if size(x.derivatives, 2) ~= size(y.derivatives, 2)
        [x,y] = binary_ext(x,y);
      end
      n_x = numel(x.values);
      n_y = numel(y.values);
      x.values = x.values - y.values;
      if n_x == n_y % numel(x) == numel(y)
        x.derivatives = x.derivatives - y.derivatives;
      elseif (n_y == 1) || (n_x == 1)
        x.derivatives = bsxfun(@minus, x.derivatives, y.derivatives);
      else
        error('Implicit expansion is not supported yet!');
      end
    else
      if numel(y) > 1 && numel(x.values) == 1
        x.derivatives = repmat(x.derivatives, numel(y), 1);
      end
      x.values = x.values - y;
    end
  else
    if numel(y) == 1 && numel(x) > 1
      y.derivatives = repmat(y.derivatives, numel(x), 1);
    end
    y.values = x - y.values ;
    y.derivatives = -y.derivatives;
    x = y;
  end
end
