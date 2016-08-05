function varargout = size(x, varargin)
% by SeHyoun Ahn, July 2016
if nargout<=1
    varargout = {size(x.values,varargin{:})};
else
    aux = size(x.values,varargin{:});
    varargout = mat2cell(aux(1:nargout),1,ones(nargout,1));
end
