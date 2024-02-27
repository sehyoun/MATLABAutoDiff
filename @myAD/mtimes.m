function x = mtimes(x, y)
% by SeHyoun Ahn, Jan 2016
% Note that in the original package by Martin Fink, this redirected to element-wise multiplication.

  if isa(x,'myAD')
    [n,m] = size(x.values);
    if isa(y,'myAD')
      if size(x.derivatives, 2) ~= size(y.derivatives, 2)
        [x,y] = binary_ext(x,y);
      end
      if max(m,n)==1
        x.derivatives = valXder(y.values(:),x.derivatives) + y.derivatives*x.values;
        x.values = x.values*y.values;
      elseif numel(y.values)==1
        x.derivatives = x.derivatives*y.values + valXder(x.values(:),y.derivatives);
        x.values = x.values*y.values;
      elseif (size(y,1)==m)
        if (size(y,2)==1)
          x.derivatives = sparse(x.values)*y.derivatives + matdrivXvecval(x.derivatives, y.values);
        else
  	      x.derivatives = matvalXmatder(x.values, y.derivatives) + matderXmatval(x.derivatives, y.values);
        end
        x.values = x.values*y.values;
      else
        error('Check that the dimensions match');
      end
    else
      if max(m,n) == 1
        x.derivatives = valXder(y(:),x.derivatives);
        x.values = x.values*y;
      elseif numel(y) == 1
        x.derivatives = x.derivatives*y;
        x.values = x.values*y;
      elseif size(y,1)==m
        if (size(y,2) == 1)
          x.derivatives = matdrivXvecval(x.derivatives, y);
        else
          x.derivatives = matderXmatval(x.derivatives, y);
        end
        x.values = x.values*y;
      else
        error('Check that the dimensions match');
      end
    end
  else
    [n,m]=size(x);
    if max(m,n)==1
      y.derivatives = y.derivatives*x;
      y.values = x*y.values;
      x = y;
    elseif numel(y.values) == 1
      y.derivatives = valXder(x(:),y.derivatives);
      y.values = x*y.values;
      x = y;
    elseif size(y,1)==m
  	  if size(y,2)==1
  	    y.derivatives = sparse(x)*y.derivatives;
      else
  	    y.derivatives = matvalXmatder(x, y.derivatives);
      end
      y.values = x*y.values;
      x = y;
    else
      error('Check that the dimensions match');
    end
  end
end
