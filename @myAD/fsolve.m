function x = fsolve(h,x0,varargin)
    tmp=@(z) h([z;getvalue(varargin{1})]);
    x.values=fsolve(tmp,x0,varargin{2:end});
    z=myAD([x.values;getvalue(varargin{1})]);
    tmp_deriv=getderivs(h(z));
    x.derivatives=-tmp_deriv(1)\tmp_deriv(2);
