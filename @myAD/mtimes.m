function x = mtimes(x, y)
% by SeHyoun Ahn, Jan 2016

% Note that in the original package by Martin Fink, this redirected to element-wise multiplication

if isa(x,'myAD')
    [n,m]=size(x.values);
    if isa(y,'myAD')
        if m>1 && size(y,1)==m
            if (size(y,2)==1)
                x.derivatives = sparse(x.values)*y.derivatives...
                                +matdrivXvecval(x.derivatives,y.values);
            else
	      x.derivatives = matvalXmatder(sparse(x.values),dertransp(y.derivatives,m)) ...
		+ dertransp(matvalXmatder(sparse(y.values'),x.derivatives),size(y,2));
            end
            x.values = x.values*y.values;
        elseif max(m,n)==1
            x.derivatives= valXder(y.values(:),x.derivatives) + y.derivatives*x.values;
            x.values=x.values*y.values;
        elseif numel(y.values)==1
            x.derivatives = x.derivatives*y.values + valXder(x.values(:),y.derivatives);
            x.values=x.values*y.values;
        else
            error('Check that the dimensions match');
        end
    else
        if m>1 && size(y,1)==m
            x.derivatives = dertransp(matvalXmatder(sparse(y'),x.derivatives),size(y,2));
            x.values = x.values*y;
        elseif max(m,n)==1
            x.derivatives= valXder(y(:),x.derivatives);
            x.values=x.values*y;
        elseif numel(y)==1
            x.derivatives=x.derivatives*y;
            x.values=x.values*y;
        else
            error('Check that the dimensions match');
        end
    end
else
    [n,m]=size(x);
    if m>1 && size(y,1)==m
		if size(y,2)==1
			y.derivatives = sparse(x)*y.derivatives;
        else
	      y.derivatives = matvalXmatder(sparse(x),dertransp(y.derivatives,m));
        end
        y.values = x*y.values;
        x=y;
    elseif max(m,n)==1
        y.derivatives=y.derivatives*x;
        y.values=x*y.values;
        x=y;
    elseif numel(y.values)==1
        y.derivatives=valXder(x(:),y.derivatives);
        y.values=x*y.values;
        x=y;
    else
        error('Check that the dimensions match');
    end
end
