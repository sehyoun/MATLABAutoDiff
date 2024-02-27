function x = vertcat(varargin)
  % by SeHyoun Ahn, Jan 2016

  i = 1;
  while (~isa(varargin{i}, 'myAD'))
    i=i+1;
  end
  y=vertcat(varargin{1:i-1});

  x = varargin{i};
  l = size(x.derivatives,2);
  if (i>1)
    aux = size(y);
    n = numel(y);
    locs = reshape(1:n,aux);
    aux = size(x.values);
    x.values = [y; x.values];
    [ii, jj, vv] = find(x.derivatives);
    ii = ii + n;
    locs = [locs; n+reshape(1:prod(aux),aux)];
    n = n+prod(aux);
  else
    aux = size(x.values);
    [ii, jj, vv] = find(x.derivatives);
    locs = reshape(1:prod(aux),aux);
    n = prod(aux);
  end

  ii_cell = cell(nargin-i+1,1);
  ii_cell{1} = ii(:);
  jj_cell = cell(nargin-i+1,1);
  jj_cell{1} = jj(:);
  vv_cell = cell(nargin-i+1,1);
  vv_cell{1} = vv(:);
  for j = i+1:nargin
    if isa(varargin{j}, 'myAD')
      x.values = [x.values; varargin{j}.values];
      [ii_curr, jj_curr, vv_curr] = find(varargin{j}.derivatives);
      % ii = [ii(:); ii_curr(:) + n];
      % jj = [jj(:); jj_curr(:)];
      % vv = [vv(:); vv_curr(:)];
      ii_cell{j-i+1} = ii_curr(:) + n;
      jj_cell{j-i+1} = jj_curr(:);
      vv_cell{j-i+1} = vv_curr(:);
      aux = size(varargin{j}.values);
    elseif (~isempty(varargin{j}))
      x.values = [x.values; varargin{j}];
      % x.derivatives = [x.derivatives;sparse(numel(varargin{j}),l)];
      aux = size(varargin{j});
    end
    locs = [locs; n+reshape(1:prod(aux),aux)];
    n = n+prod(aux);
  end
  locs = locs(:);
  ii = cell2mat(ii_cell);
  jj = cell2mat(jj_cell);
  vv = cell2mat(vv_cell);

  % x.derivatives = sparse(locs(ii), jj, vv, n, l);

  [~, inds] = sort(locs);
  x.derivatives = sparse(inds(ii), jj, vv, n, l);
  % x.derivatives = x.derivatives(locs(:),:);
end
