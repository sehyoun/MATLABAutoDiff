function h = scatter3(varargin)
% by SeHyoun Ahn, July 2020

  varargin = clean_ad(varargin{:});
  if nargout > 0
    h = scatter3(varargin{:});
  else
    scatter3(varargin{:});
  end
end
