function [output] = matvalXmatder(A, B)
  % Compute A*dB/dx
  % by SeHyoun Ahn, July 2018

  [Arow, Acol, Aval] = find(A);
  [nrow, ninter] = size(A);
  Arow = Arow(:)';
  Acol = Acol(:)';
  Aval = Aval(:)';

  [Brow, Bcol, Bval] = find(B);
  [ncol, nderiv] = size(B);
  Brow = Brow(:);
  Bcol = Bcol(:);
  Bval = Bval(:);

  ncol = ncol/ninter;

  for iter_overlap = ninter:-1:1
    ind_A = (Acol == iter_overlap);
    ind_B = (mod(Brow-1, ninter) == iter_overlap-1);

    n_inter_A = sum(ind_A);
    n_inter_B = sum(ind_B);

    row_stack{iter_overlap} = Arow(ind_A) + nrow*floor((Brow(ind_B)-1)/ninter);
    col_stack{iter_overlap} = ones(1, n_inter_A).*Bcol(ind_B);
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
