function x = subsref(x, S)
% Edited by SeHyoun Ahn, July 2016
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at

aux=size(x.values);
locs=reshape(1:prod(aux),aux);
locs=locs(S.subs{:});
x.derivatives = x.derivatives(locs(:),:);
x.values = x.values(S.subs{:});