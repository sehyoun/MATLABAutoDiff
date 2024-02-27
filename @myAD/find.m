function [i,j,v] = find(A);
% by SeHyoun Ahn, Oct 2017

% This function does not support all use cases of find.

  val = getvalues(A);
  der = getderivs(A);

  if nargout > 1
    [i,j,retval] = find(val);
    m = size(val,1);
    retder = der((j-1)*m + i, :);
    v = myAD(retval,retder);
  else
    error('find is only supported for input in sparse representation');
  end
end
