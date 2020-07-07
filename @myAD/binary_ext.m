function [x, y] = binary_ext(x, y)
  % INTERNAL FUNCTION
  % by SeHyoun Ahn, April 2020
  l_x = size(x.derivatives, 2);
  l_y = size(y.derivatives, 2);
  if l_x == l_y
    return;
  end

  warning('AutoDiff:autoext', 'Dimensionality of the derivatives are automatically extended. Working with implicitly functions result in this. Otherwise, make sure that this is the intended behavior. To the turnoff the warning, run <warning(''off'',''AutoDiff:autoext'')>.');

  if l_x > l_y
    y.derivatives(end, l_x) = 0;
  else
    x.derivatives(end, l_y) = 0;
  end
end
