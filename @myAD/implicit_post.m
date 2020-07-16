function x = implicit_post(x0, eqns)
  % by SeHyoun Ahn, April 2020
  if isa(x0, 'myAD')
    x0 = getvalues(x0);
  end

  if isa(eqns, 'myAD')
    eqns = getderivs(eqns);
  end
  n_x = size(eqns, 1);

  eqns = -eqns(:, end-n_x+1:end)\eqns(:, 1:end-n_x);
  x = myAD(x0, eqns);
end
