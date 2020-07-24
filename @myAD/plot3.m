function h = plot3(varargin)
% by SeHyoun Ahn, July 2020

  varargin = clean_ad(varargin{:});
  if nargout > 0
    h = plot3(varargin{:});
  else
    plot3(varargin{:});
  end
end
