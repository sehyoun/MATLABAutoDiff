function y = isnan(x)
  % In Package myAD - Automatic Differentiation
  % by Martin Fink, June 2006
  % martinfink 'at' gmx.at

  % Updated by SeHyoun, Jan 2019
  y = isnan(x.values) | reshape(isnan(sum(x.derivatives,2)), size(x.values));
end
