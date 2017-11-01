function varargout = max(x,varargin)
% Edited by SeHyoun Ahn, Jan 2016

% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at
if nargin > 2
    if isa(varargin{2},'double')
        direction = varargin{2};
        z=x;

        [z.values,idx] = max(x.values,varargin{:});

        aux = size(x.values);
        n_pre = prod(aux(1:direction-1));
        n_cur = prod(aux(1:direction));
        if direction > 2
            locs = bsxfun(@plus,reshape(1:n_pre,aux(1:direction-1)),(idx-1)*n_pre);
        elseif direction == 2
            locs = bsxfun(@plus,(1:n_pre)',(idx-1)*n_pre);
        else
            locs = bsxfun(@plus,(1:n_pre),(idx-1)*n_pre);
        end
        locs = bsxfun(@plus,locs,reshape(0:n_cur:prod(aux)-1,[ones(1,direction),aux(direction+1:end)]));

        z.derivatives = x.derivatives(locs(:),:);
        varargout{2} = idx;
    else
        if isa(varargin{1},'myAD')
            idx = x.values<varargin{1}.values;
            z=x;
            z.values = x.values.*(1-idx)+idx.*varargin{1}.values;
            z.derivatives(idx(:),:) = varargin{1}.derivatives(idx(:),:);
            warning('AutoDiff:maxmin','There is an ambiguity in what the derivative should be when the values are equal. This is resolved by picking the derivatives of the first one.  To the turnoff warning, run <warning(''off'',''AutoDiff:maxmin'')>.');
        else
            if isempty(varargin{1})
                z=x;
                [n,m]=size(x.values);
                [z.values,idx]=max(x.values,varargin{:});
                z.derivatives = x.derivatives(n*(0:(m-1))+idx,:);
                varargout{2} = idx;
            else
                idx = x.values<varargin{1};
                z=x;
                z.values = x.values.*(1-idx)+idx.*varargin{1};
                z.derivatives(idx(:),:) = 0;
                warning('AutoDiff:maxmin','There is an ambiguity in what the derivative should be when the values are equal. This is resolved by picking the derivatives of the first one.  To the turnoff warning, run <warning(''off'',''AutoDiff:maxmin'')>.');
            end
        end
    end
elseif nargin==2
    if isa(varargin{1},'myAD')
        idx = x.values<varargin{1}.values;
        z=x;
        z.values = x.values.*(1-idx)+idx.*varargin{1}.values;
        z.derivatives(idx(:),:) = varargin{1}.derivatives(idx(:),:);
        warning('AutoDiff:maxmin','There is an ambiguity in what the derivative should be when the values are equal. This is resolved by picking the derivatives of the first one.  To the turnoff warning, run <warning(''off'',''AutoDiff:maxmin'')>.');
    else
        idx = x.values<varargin{1};
        z=x;
        z.values = x.values.*(1-idx)+idx.*varargin{1};
        z.derivatives(idx(:),:) = 0;
        warning('AutoDiff:maxmin','There is an ambiguity in what the derivative should be when the values are equal. This is resolved by picking the derivatives of the first one.  To the turnoff warning, run <warning(''off'',''AutoDiff:maxmin'')>.');
    end
else
    z=x;
    [n,m]=size(x.values);
    [z.values,idx]=max(x.values,varargin{:});
    z.derivatives = x.derivatives(n*(0:(m-1))+idx,:);
    varargout{2} = idx;
end
varargout{1} = z;
