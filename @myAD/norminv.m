function x = norminv(x)
  % by SeHyoun Ahn, May 2022
  
  if nargin > 1
    error('General form of norminv has not been implemented');
  end
  x.values = norminv(x.values);
  x.derivatives = valXder(1./normpdf(x.values), x.derivatives);
end

