function x = subsref(x, S)
  % Edited by SeHyoun Ahn, April 2020
  % Edited by SeHyoun Ahn, July 2016
  % Edited by SeHyoun Ahn, Jan 2016

  % In Package myAD - Automatic Differentiation
  % by Martin Fink, May 2007
  % martinfink 'at' gmx.at

  if length(S) > 1
    subed = builtin('subsref', x, S(1:end-1));
    x = subsref(subed, S(end));
    return;
  end

  if strcmp(S.type,'()')
    if (length(S.subs) == 1) && strcmp(S.subs{1},':')
      x.values = x.values(:);
      return;
    end

    if isa(x.values, 'diagmat')
      x.values = sparse(x.values);
    end
    aux=size(x.values);
    if ~issparse(x.values)
      if (length(S.subs) == 1)
        x.derivatives = x.derivatives(S.subs{:},:);
        x.values = x.values(S.subs{:});
      else
        locs=reshape(1:prod(aux),aux);
        locs=locs(S.subs{:});
        locs=locs(:);

        % indexer = sparse(1:length(locs), locs, ones(length(locs),1), length(locs), size(x.derivatives, 1));
        % x.derivatives = indexer*x.derivatives;

        % [ii,jj,vv] = find(x.derivatives);
        % [ind, locator] = ismember(ii, locs);
        % jj = jj(ind);
        % vv = vv(ind);
        % x.derivatives = sparse(locator(ind), jj, vv, length(locs), size(x.derivatives,2));

        x.derivatives = x.derivatives(locs(:),:);
        x.values = x.values(S.subs{:});
      end
    else
      x.values = x.values(S.subs{:});
      if aux(2)==1
        x.derivatives = x.derivatives(S.subs{1},:);
      elseif aux(1)==1
        x.derivatives = x.derivatives(S.subs{1},:);
      else
        I = aux(1);
        if strcmp(S.subs{1},':')
          S.subs{1} = 1:aux(1);
        end
        if strcmp(S.subs{2},':')
          S.subs{2} = 1:aux(2);
        end

        if islogical(S.subs{1})
          S.subs{1} = find(S.subs{1});
        end
        if islogical(S.subs{2})
          S.subs{2} = find(S.subs{2});
        end

        S.subs{1} = mod(S.subs{1}-1,aux(1))+1;

        % func1 = @(i) find(i==S.subs{1});
        % func2 = @(i) find(i==S.subs{2});

        n = length(S.subs{1});
        m = length(S.subs{2});
        l = size(x.derivatives,2);

        func_1 = zeros(size(x,1),1);
        func_2 = zeros(size(x,2),1);

        func_1(S.subs{1}) = 1:n;
        func_2(S.subs{2}) = 1:m;

        [ii,jj,vv] = find(x.derivatives);
        locs = ismember(ceil(ii/I), S.subs{2});
        ii = ii(locs);
        jj = jj(locs);
        vv = vv(locs);


        locs = ismember(mod(int64(ii-1),int64(I))+1, S.subs{1});
        ii = ii(locs);
        jj = jj(locs);
        vv = vv(locs);
        %ii = arrayfun(func1,mod(ii(locs),I)) + (arrayfun(func2,ceil(ii(locs)/I))-1)*length(S.subs{1});
        ii = func_1(mod(int64(ii-1),int64(I))+1) + (func_2(ceil(ii/I))-1)*n;

        if length(vv) > 0 && all(diff(jj) >= 0)
          x.derivatives = sparse_bypass(ii, jj, vv, n*m, l);
        else
          x.derivatives = sparse(ii, jj, vv, n*m, l);
        end
      end
    end
  else
    x = builtin('subsref', x, S);
  end
end
