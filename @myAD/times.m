function x = times(x,y)
  % Edited by SeHyoun Ahn, Jan 2016

  % In Package myAD - Automatic Differentiation
  % by Martin Fink, May 2007
  % martinfink 'at' gmx.at
  if isa(x, 'myAD')
    if isa(y, 'myAD')
      [x,y] = binary_ext(x,y);
      if numel(y.values) == 1
        x.derivatives = y.values.*x.derivatives;
      else
        x.derivatives = valXder(y.values(:), x.derivatives);
      end
      if numel(x.values) == 1
        x.derivatives = x.derivatives + x.values(:).*y.derivatives;
      else
        x.derivatives = x.derivatives + valXder(x.values(:),y.derivatives);
      end
      x.values = x.values.*y.values;
    else
      if numel(y) == 1
        x.derivatives = y.*x.derivatives;
      else
        x.derivatives = valXder(y(:),x.derivatives);
      end
      x.values = x.values.*y;
    end
  else
    if numel(x) == 1
      y.derivatives = x.*y.derivatives;
    else
      y.derivatives = valXder(x(:) + 0.*y.values(:),y.derivatives);
    end
    y.values = x.*y.values;
    x = y;
  end
end
