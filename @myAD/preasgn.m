function y = preasgn(y, x)
  % by SeHyoun Ahn, July 2020

  n_l = size(getderivs(x), 2);
  y = myAD(y, sparse(numel(y), n_l));
end
