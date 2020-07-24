function h = plot(varargin)
% by SeHyoun Ahn, July 2020

  varargin = clean_ad(varargin{:});
  if nargout > 0
    h = plot(varargin{:});
  else
    plot(varargin{:});
  end
end
