function x = sinh(x)
% Written by SeHyoun Ahn, Jan 2016

x.derivatives = valXder(cosh(x.values), x.derivatives);
x.values = sinh(x.values);
