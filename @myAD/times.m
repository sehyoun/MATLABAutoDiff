function x = times(x,y)
  % Edited by SeHyoun Ahn, Jan 2016

  % In Package myAD - Automatic Differentiation
  % by Martin Fink, May 2007
  % martinfink 'at' gmx.at
  if isa(x, 'myAD')
    if isa(y, 'myAD')
      x.derivatives = valXder(y.values(:),x.derivatives) + valXder(x.values(:),y.derivatives);
      x.values = x.values.*y.values;
    else
      x.derivatives = valXder(y(:),x.derivatives);
      x.values = x.values.*y;
    end
  else
    y.derivatives = valXder(x(:),y.derivatives);
    y.values = x.*y.values;
    x = y;
  end
end
