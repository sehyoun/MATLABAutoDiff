function [output] = matderXmatval(A, B)
  % Compute dA/dx*B
  % by SeHyoun Ahn, July 2018

  [Arow, Acol, Aval] = find(A);
  [nrow, nderiv] = size(A);
  Arow = Arow(:)';
  Acol = Acol(:)';
  Aval = Aval(:)';

  [Brow, Bcol, Bval] = find(B);
  [ninter, ncol] = size(B);
  Brow = Brow(:);
  Bcol = Bcol(:);
  Bval = Bval(:);

  nrow = nrow / ninter;

  for iter_overlap = ninter:-1:1
    ind_A = (ceil(Arow/nrow) == iter_overlap);
    ind_B = (Brow == iter_overlap);

    n_inter_A = sum(ind_A);
    n_inter_B = sum(ind_B);

    row_stack{iter_overlap} = mod(Arow(ind_A)-1, nrow) + 1 + nrow*(Bcol(ind_B)-1);
    col_stack{iter_overlap} = ones(n_inter_B, 1).*Acol(ind_A);
    val_stack{iter_overlap} = Aval(ind_A).*Bval(ind_B);

    row_stack{iter_overlap} = row_stack{iter_overlap}(:);
    col_stack{iter_overlap} = col_stack{iter_overlap}(:);
    val_stack{iter_overlap} = val_stack{iter_overlap}(:);
  end

  row_stack = cell2mat(row_stack(:));
  col_stack = cell2mat(col_stack(:));
  val_stack = cell2mat(val_stack(:));

  output = sparse(row_stack, col_stack, val_stack, nrow*ncol, nderiv);
end
