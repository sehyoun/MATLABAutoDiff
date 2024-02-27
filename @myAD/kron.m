function output = kron(A,B)
  % by SeHyoun Ahn, May 2023

  % Can be improved for speed depending on applications

  [n_a, m_a] = size(A);
  [n_b, m_b] = size(B);

  B_ext = repmat(B, n_a, m_a);

  [ii, jj, vv] = find(A);
  ii = repmat(ii(:)-1, 1, n_b, m_b)*n_b + reshape((1:n_b), 1, n_b, 1);
  jj = repmat(jj(:)-1, 1, n_b, m_b)*m_b + reshape((1:m_b), 1, 1, m_b);
  vv = repmat(vv(:), 1, n_b*m_b);
  A_ext = sparse(ii(:), jj(:), vv(:), n_a*n_b, m_a*m_b);

  output = A_ext.*B_ext;
end
