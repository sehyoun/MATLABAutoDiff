function [U, T] = schur_helper(A, U, T)
% by SeHyoun Ahn, July 2019

  A_dot = getderivs(A);
  A = getvalues(A);

  n_row = size(A,1);
  n_k = size(U, 2);
  n_div = size(A_dot, 2);

  if max(max(abs(A*U - U*T))) > 1e-13
    warning('This is a proper subspace.');
  end

  ind_keep = find((1:n_k) - (1:n_k)' >= 0);

  A_to_flip = zeros(n_row*n_k + n_k*(n_k+1)/2, n_row*n_k + n_k*(n_k+1)/2);
  b_to_flip = zeros(n_row*n_k + n_k*(n_k+1)/2, n_div);

  T_mat = -kron(speye(n_k), U);

  A_to_flip(1:n_row*n_k,:) = [kron(speye(n_k), A) - kron(T.', speye(n_row)), T_mat(:, ind_keep)];

  full_prod = repmat(U.', n_k, 1);
  [ii, jj, vv] = find(full_prod);
  full_prod = sparse(ii, jj + mod(ii-1, n_k).*n_row, vv, n_k*n_k, n_row*n_row);
  full_prod = kron(speye(n_k), U.') + full_prod;
  % full_prod = kron(speye(n_k), U.') + kron(U.', speye(n_k));
  A_to_flip(n_row*n_k+1:end, 1:n_row*n_k) = full_prod(ind_keep, :);

  [i,j,v] = find(A_dot);
  actual_j = ceil(i/n_row);
  actual_i = mod(i-1, n_row)+1;
  v = -v.*U(actual_j, :);
  i = actual_i + n_row*(0:n_k-1);
  j = repmat(j, 1, n_k);
  b_to_flip(1:n_row*n_k, :) = sparse(i(:),j(:),v(:), n_row*n_k, n_div);

  [u_A, s_A, v_A] = svd(A_to_flip);
  big_d = sum(diag(s_A) > 1e-13);
  u_A = u_A(:, 1:big_d);
  s_A = s_A(1:big_d, 1:big_d);
  v_A = v_A(:, 1:big_d);
  result = v_A*(s_A\(u_A'*b_to_flip));

  U = myAD(U, result(1:n_row*n_k, :));
  [i,j,v] = find(result(n_row*n_k+1:end, :));
  T = myAD(T, sparse(ind_keep(i), j, v, n_k*n_k, n_div));
end

