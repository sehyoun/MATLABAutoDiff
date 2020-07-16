function [AA, BB, Q, Z] = qz_helper(A, B, AA, BB, Q, Z)
% by SeHyoun Ahn, July 2019

  A_dot = getderivs(A);
  A = getvalues(A);

  B_dot = getderivs(B);
  B = getvalues(B);

  n_row = size(A,1);
  n_k = size(Q, 2);
  n_div = size(A_dot, 2);
  n_points = 2*n_row*n_k + n_k*(n_k+1);

  Q = Q';


  ind_keep = find((1:n_k) - (1:n_k)' >= 0);

  A_to_flip = zeros(n_points, n_points);
  b_to_flip = zeros(n_points, n_div);

  T_mat = -kron(speye(n_k), Q);
  T_mat = T_mat(:, ind_keep);
  A_to_flip(1:n_row*n_k,:) = [-kron(AA.', speye(n_row)), kron(speye(n_k), A), T_mat, 0*T_mat];
  A_to_flip(n_row*n_k + 1: 2*n_row*n_k,:) = [-kron(BB.', speye(n_row)), kron(speye(n_k), B), 0*T_mat, T_mat];

  full_prod = kron(speye(n_k), Q.') + kron(Q.', speye(n_k));
  A_to_flip(2*n_row*n_k+1: 2*n_row*n_k+n_k*(n_k+1)/2, 1:n_row*n_k) = full_prod(ind_keep, :);

  full_prod = kron(speye(n_k), Z.') + kron(Z.', speye(n_k));
  A_to_flip(2*n_row*n_k+n_k*(n_k+1)/2+1:end, n_row*n_k+1:2*n_row*n_k) = full_prod(ind_keep, :);

  [i,j,v] = find(A_dot);
  actual_j = ceil(i/n_row);
  actual_i = mod(i-1, n_row)+1;
  v = -v.*Z(actual_j, :);
  i = actual_i + n_row*(0:n_k-1);
  j = repmat(j, 1, n_k);
  b_to_flip(1:n_row*n_k, :) = sparse(i(:),j(:),v(:), n_row*n_k, n_div);

  [i,j,v] = find(B_dot);
  actual_j = ceil(i/n_row);
  actual_i = mod(i-1, n_row)+1;
  v = -v.*Z(actual_j, :);
  i = actual_i + n_row*(0:n_k-1);
  j = repmat(j, 1, n_k);
  b_to_flip(n_row*n_k+1:2*n_row*n_k, :) = sparse(i(:),j(:),v(:), n_row*n_k, n_div);

  [u_A, s_A, v_A] = svd(A_to_flip);
  big_d = sum(diag(s_A) > 1e-13);
  u_A = u_A(:, 1:big_d);
  s_A = s_A(1:big_d, 1:big_d);
  v_A = v_A(:, 1:big_d);
  result = v_A*(s_A\(u_A'*b_to_flip));

  Q = myAD(Q, result(1:n_row*n_k, :));
  Q = Q';
  Z = myAD(Q, result(n_row+n_k+1:2*n_row*n_k,:));

  [i,j,v] = find(result(2*n_row*n_k+1:2*n_row*n_k+n_k*(n_k+1)/2, :));
  AA = myAD(AA, sparse(ind_keep(i), j, v, n_k*n_k, n_div));

  [i,j,v] = find(result(2*n_row*n_k+n_k*(n_k+1)/2+1:end, :));
  BB = myAD(BB, sparse(ind_keep(i), j, v, n_k*n_k, n_div));
end

