function x = erf(x)
  % by SeHyoun Ahn, July 2019

  x.derivatives = valXder(2/sqrt(pi).*exp(-x.values(:).^2), x.derivatives);
  x.values = erf(x.values);
end
