function [U, S, V] = svd_helper(A, U, S, V)
% by SeHyoun Ahn, July 2019

  A_dot = getderivs(A);
  A = getvalues(A);

  n_row = size(A,1);
  n_col = size(A,2);
  n_rank = size(U, 2);
  n_u = n_row*n_rank;
  n_v = n_col*n_rank;
  n_div = size(A_dot, 2);
  n_points = n_row*n_rank + n_col*n_rank + n_rank;

  if max(max(abs(A - U*S*V'))) > 1e-13
    warning('This is not a svd.');
  end


  ind_keep = find((1:n_rank) - (1:n_rank)' >= 0);
  A_to_flip = zeros(n_points, n_points);
  b_to_flip = zeros(n_points, n_div);
  b_to_flip(1:n_row*n_col, :) = A_dot;


  S = diag(S);
  % U_dot_S_V
  A_to_flip(1:n_row*n_col, 1:n_row*n_rank) = kron(S'.*V, speye(n_row));
  i = repmat((1:n_row*n_rank)', n_col, 1);
  j = repmat((1:n_col*n_rank), n_row, 1);
  j = j(:);
  U_tmp = U.*S';
  v = U_tmp(mod(i-1, n_row)+1  +  (ceil(j/n_col)-1)*n_row);
  % U_S_V_dot
  A_to_flip(1:n_row*n_col, n_row*n_rank + n_rank+1:end) = sparse(i, j, v, n_row*n_col, n_rank*n_col);
  % U_S_dot_V
  A_to_flip(1:n_row*n_col, n_row*n_rank + 1 :n_row*n_rank + n_rank) = kron(ones(n_col, 1), U).*kron(V, ones(n_row,1));

  full_prod_U = kron(speye(n_rank), U.') + kron(U.', speye(n_rank));
  full_prod_V = kron(speye(n_rank), V.') + kron(V.', speye(n_rank));
  A_to_flip(n_row*n_col+1: n_row*n_col + n_rank*(n_rank+1)/2, 1:n_row*n_rank)= full_prod_U(ind_keep, :);
  A_to_flip(n_row*n_col + n_rank*(n_rank+1)/2+1:end, n_row*n_rank + n_rank+1:end) = full_prod_V(ind_keep, :);

  [u_A, s_A, v_A] = svd(A_to_flip);
  big_d = sum(diag(s_A) > 1e-13);
  u_A = u_A(:, 1:big_d);
  s_A = s_A(1:big_d, 1:big_d);
  v_A = v_A(:, 1:big_d);
  result = v_A*(s_A\(u_A'*b_to_flip));

  U = myAD(U, result(1:n_row*n_rank, :));
  S = myAD(S, result(n_row*n_rank + 1: n_row*n_rank+n_rank, :));
  S = spdiags(S, 0, n_rank, n_rank);
  V = myAD(V, result(n_row*n_rank+n_rank+1:end, :));
end

