function [output] = matdrivXvecval(A,b)
  % Compute the dA/dt * b
  % by SeHyoun Ahn, July 2018

  [row,col,val] = find(A);
  [nrow,nderiv] = size(A);
  m = length(b);
  nrow = nrow/m;
  val = val.*b(ceil(row/nrow));
  output = sparse(mod(row-1, nrow)+1, col, val, nrow, nderiv);
end
