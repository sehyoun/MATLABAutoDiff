function x=spdiags(x,varargin)
    x.values=spdiags(x.values,varargin{1},varargin{2},varargin{3});
    l=size(x.derivatives,2);
    tmp=sparse(varargin{2}*varargin{3},l);
    if varargin{1}<0
        tmp(1-varargin{1}:varargin{2}+1:varargin{2}*(varargin{3}+varargin{1}),:)=x.derivatives(1:varargin{2}+varargin{1},:);
    else
        tmp(varargin{1}*varargin{2}+1:varargin{2}+1:varargin{2}*varargin{3}-varargin{1},:)=x.derivatives(varargin{1}+1:varargin{2},:);
    end
    x.derivatives=tmp;