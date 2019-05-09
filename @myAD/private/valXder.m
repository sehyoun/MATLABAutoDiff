function der = valXder(val, der)
  % Computes different instances of v*dx
  %

  val = val(:);
  n_val = length(val);
  [nrow, ncol] = size(der);

  [row, col, derval] = find(der);
  row = row(:);
  col = col(:);
  derval = derval(:);


  if (nrow == 1 & n_val > 1)
    row = repmat((1:n_val)', length(col), 1);
    col = repmat(col', n_val, 1);
    val = val(:).*derval';
    der = sparse(row, col(:), val(:), n_val, ncol);
  elseif n_val == 1
    der = sparse(row, col, derval.*val, nrow, ncol);
  else
    der = sparse(row, col, derval.*val(row), nrow, ncol);
  end
end
