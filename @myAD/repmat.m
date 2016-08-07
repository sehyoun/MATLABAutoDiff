function x=repmat(x,varargin)
% by SeHyoun Ahn, Jan 2016

if ~issparse(x.values)
    aux=size(x.values);
    p=repmat(reshape(1:prod(aux),aux),varargin{:});
    x.derivatives=x.derivatives(p(:),:);
    x.values=repmat(x.values,varargin{:});
else
    [i,j,v] = find(x.derivatives);
    l = size(x.derivatives,2);
    [n,m] = size(x.values);
    var1 = varargin{1};
    var2 = varargin{2};
    
    ix = mod(i,n);
    ix(ix==0) = n;
    iy = ceil(i/n);
    
    aux = @(ix,iy) reshape(bsxfun(@plus,ix+n*(0:var1-1)'+(iy-1)*var1*n,n*m*var1*(0:var2-1)),var1*var2,1);
    
    i = cell2mat(arrayfun(aux,ix,iy,'UniformOutput',false));
    j = kron(j,ones(var1*var2,1));
    v = kron(v,ones(var1,var2,1));
    
    x.derivatives = sparse(i,j,v,var1*var2*n*m,l);
    x.values = repmat(x.values,var1,var2);
end