function x = diag(x)
  % Edited by SeHyoun Ahn, July 2016
  % Edited by SeHyoun Ahn, Jan 2016

  % In Package myAD - Automatic Differentiation
  % by Martin Fink, May 2007
  % martinfink 'at' gmx.at

  [n_1, n_2] = size(x.values);
  n_size = min(n_1, n_2);

  ii = (1:n_size)';
  deriv_loc = ii + (ii-1).*n_1;

  x.values = diag(x.values);
  x.derivatives = x.derivatives(deriv_loc, :);
end
