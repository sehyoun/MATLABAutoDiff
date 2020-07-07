function y = subsasgn(y, S, x)
  % Edited by SeHyoun Ahn, July 2016
  % Edited by SeHyoun Ahn, Jan 2016

  % In Package myAD - Automatic Differentiation
  % by Martin Fink, May 2007
  % martinfink 'at' gmx.at

  if isa(x, 'myAD')
    if isa(y, 'myAD')
      y.values(S.subs{:}) = x.values;
      aux=size(y.values);
      locs=reshape(1:prod(aux),aux);
      locs=locs(S.subs{:});

      % n_y = size(y.derivatives,1);
      % n_loc = length(locs);
      % mid = ones(n_y,1);
      % mid(locs(:)) = 0;
      % indexer = spdiags(mid, 0, n_y, n_y);
      % y.derivatives = indexer*y.derivatives;
      % indexer = sparse(locs(:), 1:n_loc, ones(n_loc,1), n_y, n_loc);
      % y.derivatives = y.derivatives + indexer*x.derivatives;

      % [ii,jj,vv] = find(y.derivatives);
      % ind = ~ismember(ii, locs);
      % ii = ii(ind);
      % jj = jj(ind);
      % vv = vv(ind);
      % [iix,jjx,vvx] = find(x.derivatives);
      % y.derivatives = sparse([ii; locs(iix)], [jj; jjx], [vv; vvx], size(y.derivatives,1), size(y.derivatives,2));

      y.derivatives(locs(:),:) = x.derivatives;
    else
      aux=size(y);
      locs=reshape(1:prod(aux),aux);
      locs=locs(S.subs{:});
      l = size(x.derivatives,2);
      y = myAD(y,sparse(prod(aux),l));
      y.values(S.subs{:}) = x.values;
      y.derivatives(locs(:),:) = x.derivatives;
    end
  else
    aux=size(y.values);
    locs=reshape(1:prod(aux),aux);
    locs=locs(S.subs{:});
    y.values(S.subs{:}) = x;

    % n_y = size(y.derivatives,1);
    % n_loc = length(locs);
    % mid = ones(n_y,1);
    % mid(locs(:)) = 0;
    % indexer = spdiags(mid, 0, n_y, n_y);
    % y.derivatives = indexer*y.derivatives;

    % [ii,jj,vv] = find(y.derivatives);
    % ind = ~ismember(ii, locs);
    % ii = ii(ind);
    % jj = jj(ind);
    % vv = vv(ind);
    % y.derivatives = sparse(ii, jj, vv, size(y.derivatives,1), size(y.derivatives,2));

    y.derivatives(locs(:),:) = 0;
  end
end
