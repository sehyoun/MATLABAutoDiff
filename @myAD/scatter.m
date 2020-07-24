function h = scatter(varargin)
% by SeHyoun Ahn, July 2020

  varargin = clean_ad(varargin{:});
  if nargout > 0
    h = scatter(varargin{:});
  else
    scatter(varargin{:});
  end
end
